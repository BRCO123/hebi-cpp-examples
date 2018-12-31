#include <atomic>
#include <iostream>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <set>

#include "input/input_manager_mobile_io.hpp"
#include "robot/quadruped_parameters.hpp"
#include "robot/quadruped.hpp"

using namespace hebi;
using namespace Eigen;

/* state machine related needed variables */
// state definitions
enum ctrl_state_type {
  QUAD_CTRL_STAND_UP1,
  QUAD_CTRL_STAND_UP2,

  QUAD_DYNAMIC_WALK,

  CTRL_STATES_COUNT
};

// state transition variables
double startup_seconds = 1.9;
double re_pose_seconds = 3;

int main(int argc, char** argv)
{
  // INIT VARS
  bool is_quiet{}; // test where some steps go run or not

  // INIT STEP 1: init parameters (it is empty now, not used)
  QuadrupedParameters params;
  params.resetToDefaults();

  // INIT STEP 2: init input
  std::unique_ptr<input::InputManager> input(new input::InputManagerMobileIO());
  while (!input->isConnected())
  {
    static_cast<input::InputManagerMobileIO*>(input.get())->reset();
  }

  std::cout << "Found input joystick -- starting control program.\n";
  // INIT STEP 3: init robot planner
  std::unique_ptr<Quadruped> quadruped = Quadruped::create(params);
  quadruped -> setGains();

  // INIT STEP FINAL: start control state machine
  Eigen::Vector3f translation_velocity_cmd;
  translation_velocity_cmd.setZero();
  Eigen::Vector3f rotation_velocity_cmd;
  rotation_velocity_cmd.setZero();   


  // START CONTROL THREAD: this is so cool
  auto start_time = std::chrono::steady_clock::now();
  // some time for state use
  auto state_enter_time = std::chrono::steady_clock::now(); 
  auto state_curr_time = std::chrono::steady_clock::now(); 
  std::chrono::duration<double> state_run_time = std::chrono::seconds(0);
  long interval_ms = 5.0; // in milliseconds; e.g., 5 ms => 1000/5 = 200 Hz
  // http://stackoverflow.com/questions/30425772/c-11-calling-a-c-function-periodically
  std::atomic<bool> control_execute;
  control_execute.store(true, std::memory_order_release);
  int move_leg_id;

  // the main control thread
  ctrl_state_type cur_ctrl_state = QUAD_CTRL_STAND_UP1;

  auto prev_time = std::chrono::steady_clock::now();
  // Get dt (in seconds)
  std::chrono::duration<double> dt;

  // variables used in normal run test
  Eigen::Vector3d grav_vec;
  Eigen::VectorXd leg_angles;
  bool isFinished;

  // some variables used in state passive_orient
  Eigen::Matrix3d bodyR, balance_body_R, diff_body_R;
  Eigen::Vector3d euler;
  Eigen::Quaterniond q_error;
  Eigen::Matrix3d control_R = Eigen::Matrix3d::Identity();
  Eigen::AngleAxisd tmp_aa;
  double new_angle;
  Eigen::Vector3d axis_aa;
  // used in rePose and extend
  double tgt_x, tgt_y;
  bool first_time_enter = false;

  while (control_execute.load(std::memory_order_acquire)) // I have no idea why this not working after remove Q application
  // while(true)
  {    
    // Wait to make sure the loop frequency
    auto now_time = std::chrono::steady_clock::now();
    auto need_to_wait = std::max(0, (int)std::chrono::duration_cast<std::chrono::milliseconds>(prev_time + std::chrono::milliseconds(interval_ms) - now_time).count());
    std::this_thread::sleep_for(std::chrono::milliseconds(need_to_wait));

    // Get dt (in seconds)
    now_time = std::chrono::steady_clock::now();
    dt = std::chrono::duration_cast<std::chrono::duration<double>>(now_time - prev_time);
    prev_time = now_time;

    std::chrono::duration<double> elapsed_time(now_time - start_time);

    // Get joystick update, and update any relevant variables.
    input->update();
    if (input->getQuitButtonPushed())
    {
      std::cout << "program exit" << std::endl;
      std::terminate();
    }
    translation_velocity_cmd = input->getTranslationVelocityCmd();
    rotation_velocity_cmd = input->getRotationVelocityCmd();
    Eigen::Matrix3d target_body_R;

    switch (cur_ctrl_state)
    {
      case QUAD_CTRL_STAND_UP1:
      {
        state_curr_time = std::chrono::steady_clock::now();
        state_run_time = std::chrono::duration_cast<std::chrono::duration<double>>(state_curr_time - state_enter_time);

        isFinished = quadruped -> spreadAllLegs();
        if (state_run_time.count() >= startup_seconds)
        {
          cur_ctrl_state = QUAD_CTRL_STAND_UP2;
          state_enter_time = std::chrono::steady_clock::now(); 
        }
        break;
      }
      case QUAD_CTRL_STAND_UP2:
      {
        state_curr_time = std::chrono::steady_clock::now();
        state_run_time = std::chrono::duration_cast<std::chrono::duration<double>>(state_curr_time - state_enter_time);

        isFinished = quadruped -> pushAllLegs(state_run_time.count(), startup_seconds);
        if (state_run_time.count() >= startup_seconds)
        {
          quadruped -> startBodyRUpdate();
          cur_ctrl_state = QUAD_DYNAMIC_WALK;
          first_time_enter = true;
          state_enter_time = std::chrono::steady_clock::now(); 
        }
        break;
      }

      case QUAD_DYNAMIC_WALK:
      {
        state_curr_time = std::chrono::steady_clock::now();
        state_run_time = std::chrono::duration_cast<std::chrono::duration<double>>(state_curr_time - state_enter_time);
        
        // make sure the time is strictly within totalTime
        if (state_run_time.count() >= quadruped -> getTotalTime()){
          state_enter_time = std::chrono::steady_clock::now();
          first_time_enter = true;
          break;
        }

        if(first_time_enter){
          quadruped -> planDynamicGait();
          first_time_enter = false;
        }

        quadruped -> followDynamicGait(state_run_time.count());

        break;
      }

      case CTRL_STATES_COUNT:
      default:
        quadruped -> freeze();
        break;

    }
  }
  
  // control_execute.store(false, std::memory_order_release);
  return 0;
}