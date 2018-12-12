#pragma once

#include <Eigen/Dense>

class MadgwickAHRS
{
  public:
    MadgwickAHRS();
    void MadgwickAHRSupdateIMU(double gx, double gy, double gz, double ax, double ay, double az, double dt);
    void setSampleFreq(double _sampleFreq);
    
    Eigen::Quaterniond getOrientation();

  private:
    double invSqrt(double x);

    double sampleFreq;
    double q0, q1, q2, q3;   // quaternion of sensor frame relative to auxiliary frame
    double beta;             // algorithm gain

};