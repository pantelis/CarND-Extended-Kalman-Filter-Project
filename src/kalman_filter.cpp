#include <iostream>
#include <math.h>
#include "kalman_filter.h"
#include "tools.h"



using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {

    // predict the state
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::UpdateCommon(const VectorXd& y){

    const MatrixXd PHt = P_ * H_.transpose();
    const MatrixXd S = H_ * PHt + R_;
    const MatrixXd K = PHt * S.inverse();

    x_ += K * y;
    P_ -= K * H_ * P_;
}

void KalmanFilter::Update(const VectorXd &z) {

    // update the state by using Kalman Filter equations
    VectorXd z_pred_ = H_ * x_;
    VectorXd y = z - z_pred_;

    UpdateCommon(y);
}

inline void wrapAngletoPi( double& angle )
{
    double twoPi = 2.0 * M_PI;
    angle = angle - twoPi * floor( (angle+M_PI) / twoPi );
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    //update the state by using Extended Kalman Filter equations

    long x_size_ = x_.size();

    //recover predicted state parameters
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    // Map predicted state into measurement space via the measurement function
    // Convert to polar coordinates of the measurement space
    // from cartesian coordinates.

    double rho = sqrt(px*px+py*py);
    double phi = atan2(py, px);
    double rho_dot = (px*vx+py*vy)/rho;


    VectorXd z_pred_(3);
    z_pred_ << rho, phi, rho_dot;

    VectorXd y = z - z_pred_;

/*
    Imagine you're sitting in the car and looking straightforward. Your phi angle is 0.
    If you turn your head to the left 90 degrees, your phi is 90 degrees (or pi/2 rad).
    If you turn your head to the right 90 degrees, you phi is -90 degrees (- pi/2 rad).
    If you keep turning your head to the right and look back, the phi will be -180 degrees (- pi rad).
    If you turn your head to the left until you look back, the phi will be growing until 180 degrees (pi rad).

    So your phi is in the range of [-pi; pi].
   The angle 1.5*pi doesn't make sense. You know that if you keep turning your head to the left,
   the angle will be changing like 0 -> 1/8*pi ->1/2*pi -> 7/8*pi -> pi -> -7/8* pi -> -pi/2 -> -1/8*pi -> 0 -> 1/8*pi.

    From the math in kalman filter, where you sum up angles or substract angles, you might get large angle value.
    Angle normalization brings this value back in the range [-pi; pi].
*/
    wrapAngletoPi(y(1));

    UpdateCommon(y);

}
