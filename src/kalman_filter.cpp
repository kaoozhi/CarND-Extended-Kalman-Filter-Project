#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  /**
   * TODO: predict the state
   */
  // Predict state using motion model
  x_ = F_ * x_;
  // Caculate covariance matrix
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  // Caculate radar measurement matrix h'(x) (polar/cartesian conversion)
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float c1 = px*px+py*py;
  float c2 = px*vx+py*vy;
  VectorXd Hprime_(3);
  float phi = atan2(py,px);

  // Avoid division by zero
  if (fabs(c1) < 0.0001) 
  {
    c1 =0.0001;
  }
  Hprime_ << sqrt(c1),
          phi,
          c2/sqrt(c1);

  // Caculate estimation error using measurement matrix h'(x)
  VectorXd y = z - Hprime_;
  double y_phi = y(1);

  // Wrap estimation error of phi within (-pi, pi]
  if (y_phi > M_PI){y_phi=y_phi - 2*M_PI;}
  else if (y_phi <= -M_PI){y_phi=y_phi + 2*M_PI;}
  y(1) = y_phi;

  // Caculate EKF's gain using Jacobien matrix Hj
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // Update state and covariance matrix
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
