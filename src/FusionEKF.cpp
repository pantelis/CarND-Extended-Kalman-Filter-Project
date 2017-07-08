#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;


    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    Hj_ << 1, 1, 0, 0,
            1, 1, 0, 0,
            1, 1, 1, 1;


    // A-priori optimized parameters: Use noise_ax = 9 and noise_ay = 9 for the Q matrix.
    noise_ax = 9.0;
    noise_ay = 9.0;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {

        //  Initialize the state ekf_.x_ with the first measurement.

        // first measurement
        cout << "EKF: " << endl;
        // In the case of the Lidar the state will be [x, y, vx, vy]
        // In the case of the Radar the state will be [rho, phi, and Doppler]
        ekf_.x_ = VectorXd(4);

        //state covariance matrix P
        ekf_.P_ = MatrixXd(4, 4);

        ekf_.P_ << 1., 0., 0., 0.,
                0., 1., 0., 0.,
                0., 0., 100., 0.,
                0., 0., 0., 1000.;

        //the initial transition matrix F_
        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ << 1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;

        // Process noise covariance Matrix Q_
        ekf_.Q_ = MatrixXd(4, 4);
        ekf_.Q_ << 1, 0, 1, 0,
                0, 1, 0, 1,
                1, 0, 1, 0,
                0, 1, 0, 1;

        //#R(for radar) meas_rho meas_phi meas_rho_dot timestamp gt_px gt_py gt_vx gt_vy
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

            // Convert radar from polar to cartesian coordinates and initialize state.
            // x = rho * cos(phi), y = rho * sin(phi)
            ekf_.x_ << measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]),
                    measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]),
                    measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]),
                    measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);

            ekf_.H_ = Hj_;
            ekf_.R_ = R_radar_;

        }// #L(for laser) meas_px meas_py timestamp gt_px gt_py gt_vx gt_vy
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
                    0.0, 0.0;

            ekf_.H_ = H_laser_;
            ekf_.R_ = R_laser_;
        }

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    //compute the time elapsed between the current and previous measurements
    //dt is expressed in seconds
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    // Update the state transition matrix F according to the new elapsed time.
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // Update the process covariance matrix Q according to the new elapsed time
    double term4 = pow(dt, 4.) / 4.;
    double term3 = pow(dt, 3.) / 2.;
    double term2 = pow(dt, 2.);

    ekf_.Q_(0, 0) = term4 * noise_ax;
    ekf_.Q_(0, 2) = term3 * noise_ax;
    ekf_.Q_(1, 1) = term4 * noise_ay;
    ekf_.Q_(1, 3) = term3 * noise_ax;
    ekf_.Q_(2, 0) = term3 * noise_ax;
    ekf_.Q_(2, 2) = term2 * noise_ax;
    ekf_.Q_(3, 1) = term3 * noise_ay;
    ekf_.Q_(3, 3) = term2 * noise_ay;

    ekf_.Predict();
    //cout << "Finished Prediction Step" << endl;

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    //#R(for radar) meas_rho meas_phi meas_rho_dot timestamp gt_px gt_py gt_vx gt_vy
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

        // Calculate the Jacobian at the predicted state
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

        // Update the state and covariance matrices.
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);

        //cout << "Finished Radar Update Step" << endl;
    }// #L(for laser) meas_px meas_py timestamp gt_px gt_py gt_vx gt_vy
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

        // Update the state and covariance matrices.
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);

        //cout << "Finished Lidar Update Step" << endl;
    }


    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
