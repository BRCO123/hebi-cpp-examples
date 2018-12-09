#include <atomic>
#include <iostream>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <set>

#include <QtWidgets/QApplication>

#include "input/input_manager_mobile_io.hpp"
#include "robot/quadruped_parameters.hpp"
#include "robot/quadruped.hpp"

using namespace hebi;
using namespace Eigen;

/* state machine related needed variables */
// state definitions
enum ctrl_state_type {
  HEXA_CTRL_STAND_UP_PLAN,
  HEXA_CTRL_STAND_UP,
  QUAD_CTRL_NORMAL,
  CTRL_STATES_COUNT
};
// state transition variables
double startup_seconds = 4.5;

int main(int argc, char** argv)
{
  // hexapod initially start some Qt object, it seems it is just used as
  QApplication app(argc, argv);

  // INIT VARS
  bool is_quiet{}; // test where some steps go run or not

  // INIT STEP 1: init parameters (use default parameters)
  QuadrupedParameters params;
  params.resetToDefaults();

  // INIT STEP 2: init input
  std::unique_ptr<input::InputManager> input(new input::InputManagerMobileIO());
  if (!is_quiet && !input->isConnected())
  {
      std::cout << "Could not find input joystick." << std::endl;
      return 1;
  }
  // ------------ Retry a "reset" multiple times! Wait for this in a loop.
  while (is_quiet && !input->isConnected())
  {
      static_cast<input::InputManagerMobileIO*>(input.get())->reset();
  }

  std::cout << "Found input joystick -- starting control program.\n";
  // INIT STEP 3: init robot planner
  std::unique_ptr<Quadruped> quadruped = Quadruped::create(params);

  // INIT STEP FINAL: start control state machine
  // input command from joystick (hebi's input manager use vector3f, i think use vector3d would be better)
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
  std::thread control_thread([&]()
  {
    // the main control thread
    ctrl_state_type cur_ctrl_state = HEXA_CTRL_STAND_UP_PLAN;
    auto prev_time = std::chrono::steady_clock::now();
    // Get dt (in seconds)
    std::chrono::duration<double> dt = std::chrono::seconds(0);
    while (control_execute.load(std::memory_order_acquire))
    {    
      // variables init
      Eigen::Vector3d grav_vec;
      Eigen::VectorXd leg_angles;
      // Wait!
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
        app.exit();
      }
      translation_velocity_cmd = input->getTranslationVelocityCmd();
      rotation_velocity_cmd = input->getRotationVelocityCmd();

      // control state machine (do not put actual execution logic here)
      std::cout << "|Time: " << elapsed_time.count() <<  "| my current state is: " << cur_ctrl_state <<std::endl;
      switch (cur_ctrl_state)
      {
  
        case HEXA_CTRL_STAND_UP_PLAN:
          // plan a stand up trajectory
          quadruped -> planStandUpTraj(startup_seconds);
          // state transition
          cur_ctrl_state = HEXA_CTRL_STAND_UP;
          state_enter_time = std::chrono::steady_clock::now(); 
          break;

        case HEXA_CTRL_STAND_UP:
          // start up logic 
          state_curr_time = std::chrono::steady_clock::now();
          state_run_time = std::chrono::duration_cast<std::chrono::duration<double>>(state_curr_time - state_enter_time);
          quadruped -> execStandUpTraj(state_run_time.count());
          std::cout << "state: " << state_run_time.count() << std::endl;

          if (state_run_time.count() >= startup_seconds)
          {
            cur_ctrl_state = QUAD_CTRL_NORMAL;
          }
          break;

        case QUAD_CTRL_NORMAL:
          // normal state logic 
          // some debug print to show so far all implementation is correct
          // grav_vec = quadruped -> getGravityDirection();
          // std::cout << "Gravity Direction: " << grav_vec(0) << " "
          //                                    << grav_vec(1) << " "
          //                                    << grav_vec(2) << std::endl;
          // leg_angles = quadruped -> getLegJointAngles(1);
          // std::cout << "Leg angle: " << leg_angles(0) << " "
          //                            << leg_angles(1) << " "
          //                            << leg_angles(2) << std::endl;
          
          break;

        case CTRL_STATES_COUNT:
        default:
          break;

      }
      
    }
  });

  bool res = app.exec();
  control_execute.store(false, std::memory_order_release);
  control_thread.join();
  return res;
}