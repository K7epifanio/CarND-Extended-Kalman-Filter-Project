#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initialize matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix for laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix for radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  // Measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;


  Hj_ <<
    1, 1, 0, 0,
    1, 1, 0, 0,
    1, 1, 1, 1;

  ekf_.F_ = MatrixXd(4, 4);
  //initialize state transition matrix
  ekf_.F_ <<
    1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4, 4);
  //state covariance matrix 
  ekf_.P_ <<
    1, 0, 0,    0,
    0, 1, 0,    0,
    0, 0, 500, 0,
    0, 0, 0, 500;

  // Measurement matrix for radar is calculated dynamically

  

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    

    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.Q_ = MatrixXd(4,4);


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float delta_x = rho * cos(phi); //using trigonometric cos function : cos = A/H
      float delta_y = rho * sin(phi); ////using trigonometric cos function : cos = O/H
      ekf_.x_ << delta_x, delta_y, 0, 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);
      float vx = 0;
      float vy = 0;

      // Our initial x is just the measurement with velocity 0.
      ekf_.x_ << px, py, vx, vy ;


    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // time difference should be expressed in seconds hence / by 1000000
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  float noise_ax_ = 9; //acceleration which is taken as noise in this project
  float noise_ay_ = 9;  

  // update state transition matrix according to new elapsed time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Update the process noise covariance matrix Q
  ekf_.Q_ <<
    dt_4/4*noise_ax_, 0,                dt_3/2*noise_ax_, 0,
    0,                dt_4/4*noise_ay_, 0,                dt_3/2*noise_ay_,
    dt_3/2*noise_ax_, 0,                dt_2*noise_ax_,   0,
    0,                dt_3/2*noise_ay_, 0,                dt_2*noise_ay_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_); //calling the jacobian matrix
    ekf_.H_ = Hj_; //
    ekf_.R_=R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);


  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
