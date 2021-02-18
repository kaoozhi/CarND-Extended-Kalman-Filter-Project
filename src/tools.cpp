#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * Calculate the RMSE
   */
  VectorXd rmse(4);
  rmse << 0.0,0.0,0.0,0.0;
  for (int i=0; i<estimations.size(); ++i){
     VectorXd error;
     error = estimations[i] - ground_truth[i];
     
     error  = error.array()*error.array();

     rmse += error;
  }

    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * Calculate Jacobian matrix for given state.
   */
  MatrixXd Hj(3,4);
  // recover state
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;


  // check division by zero
  if (fabs(c1) < 0.0001) {

    c1 = 0.0001;
    //  Hj << 0, 0, 0, 0,
    //   0, 0, 0, 0,
    //   0, 0, 0, 0;
  }

  float c2 = sqrt(c1);
  float c3 = (c1*c2);
  // compute the Jacobian matrix
  // else {
    Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  // }
  return Hj;
}
