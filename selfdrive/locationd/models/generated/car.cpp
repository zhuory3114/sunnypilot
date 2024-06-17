#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6221613173377875958) {
   out_6221613173377875958[0] = delta_x[0] + nom_x[0];
   out_6221613173377875958[1] = delta_x[1] + nom_x[1];
   out_6221613173377875958[2] = delta_x[2] + nom_x[2];
   out_6221613173377875958[3] = delta_x[3] + nom_x[3];
   out_6221613173377875958[4] = delta_x[4] + nom_x[4];
   out_6221613173377875958[5] = delta_x[5] + nom_x[5];
   out_6221613173377875958[6] = delta_x[6] + nom_x[6];
   out_6221613173377875958[7] = delta_x[7] + nom_x[7];
   out_6221613173377875958[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1873114568696790253) {
   out_1873114568696790253[0] = -nom_x[0] + true_x[0];
   out_1873114568696790253[1] = -nom_x[1] + true_x[1];
   out_1873114568696790253[2] = -nom_x[2] + true_x[2];
   out_1873114568696790253[3] = -nom_x[3] + true_x[3];
   out_1873114568696790253[4] = -nom_x[4] + true_x[4];
   out_1873114568696790253[5] = -nom_x[5] + true_x[5];
   out_1873114568696790253[6] = -nom_x[6] + true_x[6];
   out_1873114568696790253[7] = -nom_x[7] + true_x[7];
   out_1873114568696790253[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_9120391434276283067) {
   out_9120391434276283067[0] = 1.0;
   out_9120391434276283067[1] = 0;
   out_9120391434276283067[2] = 0;
   out_9120391434276283067[3] = 0;
   out_9120391434276283067[4] = 0;
   out_9120391434276283067[5] = 0;
   out_9120391434276283067[6] = 0;
   out_9120391434276283067[7] = 0;
   out_9120391434276283067[8] = 0;
   out_9120391434276283067[9] = 0;
   out_9120391434276283067[10] = 1.0;
   out_9120391434276283067[11] = 0;
   out_9120391434276283067[12] = 0;
   out_9120391434276283067[13] = 0;
   out_9120391434276283067[14] = 0;
   out_9120391434276283067[15] = 0;
   out_9120391434276283067[16] = 0;
   out_9120391434276283067[17] = 0;
   out_9120391434276283067[18] = 0;
   out_9120391434276283067[19] = 0;
   out_9120391434276283067[20] = 1.0;
   out_9120391434276283067[21] = 0;
   out_9120391434276283067[22] = 0;
   out_9120391434276283067[23] = 0;
   out_9120391434276283067[24] = 0;
   out_9120391434276283067[25] = 0;
   out_9120391434276283067[26] = 0;
   out_9120391434276283067[27] = 0;
   out_9120391434276283067[28] = 0;
   out_9120391434276283067[29] = 0;
   out_9120391434276283067[30] = 1.0;
   out_9120391434276283067[31] = 0;
   out_9120391434276283067[32] = 0;
   out_9120391434276283067[33] = 0;
   out_9120391434276283067[34] = 0;
   out_9120391434276283067[35] = 0;
   out_9120391434276283067[36] = 0;
   out_9120391434276283067[37] = 0;
   out_9120391434276283067[38] = 0;
   out_9120391434276283067[39] = 0;
   out_9120391434276283067[40] = 1.0;
   out_9120391434276283067[41] = 0;
   out_9120391434276283067[42] = 0;
   out_9120391434276283067[43] = 0;
   out_9120391434276283067[44] = 0;
   out_9120391434276283067[45] = 0;
   out_9120391434276283067[46] = 0;
   out_9120391434276283067[47] = 0;
   out_9120391434276283067[48] = 0;
   out_9120391434276283067[49] = 0;
   out_9120391434276283067[50] = 1.0;
   out_9120391434276283067[51] = 0;
   out_9120391434276283067[52] = 0;
   out_9120391434276283067[53] = 0;
   out_9120391434276283067[54] = 0;
   out_9120391434276283067[55] = 0;
   out_9120391434276283067[56] = 0;
   out_9120391434276283067[57] = 0;
   out_9120391434276283067[58] = 0;
   out_9120391434276283067[59] = 0;
   out_9120391434276283067[60] = 1.0;
   out_9120391434276283067[61] = 0;
   out_9120391434276283067[62] = 0;
   out_9120391434276283067[63] = 0;
   out_9120391434276283067[64] = 0;
   out_9120391434276283067[65] = 0;
   out_9120391434276283067[66] = 0;
   out_9120391434276283067[67] = 0;
   out_9120391434276283067[68] = 0;
   out_9120391434276283067[69] = 0;
   out_9120391434276283067[70] = 1.0;
   out_9120391434276283067[71] = 0;
   out_9120391434276283067[72] = 0;
   out_9120391434276283067[73] = 0;
   out_9120391434276283067[74] = 0;
   out_9120391434276283067[75] = 0;
   out_9120391434276283067[76] = 0;
   out_9120391434276283067[77] = 0;
   out_9120391434276283067[78] = 0;
   out_9120391434276283067[79] = 0;
   out_9120391434276283067[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8908813278594628189) {
   out_8908813278594628189[0] = state[0];
   out_8908813278594628189[1] = state[1];
   out_8908813278594628189[2] = state[2];
   out_8908813278594628189[3] = state[3];
   out_8908813278594628189[4] = state[4];
   out_8908813278594628189[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8908813278594628189[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8908813278594628189[7] = state[7];
   out_8908813278594628189[8] = state[8];
}
void F_fun(double *state, double dt, double *out_285691373668569915) {
   out_285691373668569915[0] = 1;
   out_285691373668569915[1] = 0;
   out_285691373668569915[2] = 0;
   out_285691373668569915[3] = 0;
   out_285691373668569915[4] = 0;
   out_285691373668569915[5] = 0;
   out_285691373668569915[6] = 0;
   out_285691373668569915[7] = 0;
   out_285691373668569915[8] = 0;
   out_285691373668569915[9] = 0;
   out_285691373668569915[10] = 1;
   out_285691373668569915[11] = 0;
   out_285691373668569915[12] = 0;
   out_285691373668569915[13] = 0;
   out_285691373668569915[14] = 0;
   out_285691373668569915[15] = 0;
   out_285691373668569915[16] = 0;
   out_285691373668569915[17] = 0;
   out_285691373668569915[18] = 0;
   out_285691373668569915[19] = 0;
   out_285691373668569915[20] = 1;
   out_285691373668569915[21] = 0;
   out_285691373668569915[22] = 0;
   out_285691373668569915[23] = 0;
   out_285691373668569915[24] = 0;
   out_285691373668569915[25] = 0;
   out_285691373668569915[26] = 0;
   out_285691373668569915[27] = 0;
   out_285691373668569915[28] = 0;
   out_285691373668569915[29] = 0;
   out_285691373668569915[30] = 1;
   out_285691373668569915[31] = 0;
   out_285691373668569915[32] = 0;
   out_285691373668569915[33] = 0;
   out_285691373668569915[34] = 0;
   out_285691373668569915[35] = 0;
   out_285691373668569915[36] = 0;
   out_285691373668569915[37] = 0;
   out_285691373668569915[38] = 0;
   out_285691373668569915[39] = 0;
   out_285691373668569915[40] = 1;
   out_285691373668569915[41] = 0;
   out_285691373668569915[42] = 0;
   out_285691373668569915[43] = 0;
   out_285691373668569915[44] = 0;
   out_285691373668569915[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_285691373668569915[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_285691373668569915[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_285691373668569915[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_285691373668569915[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_285691373668569915[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_285691373668569915[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_285691373668569915[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_285691373668569915[53] = -9.8000000000000007*dt;
   out_285691373668569915[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_285691373668569915[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_285691373668569915[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_285691373668569915[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_285691373668569915[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_285691373668569915[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_285691373668569915[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_285691373668569915[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_285691373668569915[62] = 0;
   out_285691373668569915[63] = 0;
   out_285691373668569915[64] = 0;
   out_285691373668569915[65] = 0;
   out_285691373668569915[66] = 0;
   out_285691373668569915[67] = 0;
   out_285691373668569915[68] = 0;
   out_285691373668569915[69] = 0;
   out_285691373668569915[70] = 1;
   out_285691373668569915[71] = 0;
   out_285691373668569915[72] = 0;
   out_285691373668569915[73] = 0;
   out_285691373668569915[74] = 0;
   out_285691373668569915[75] = 0;
   out_285691373668569915[76] = 0;
   out_285691373668569915[77] = 0;
   out_285691373668569915[78] = 0;
   out_285691373668569915[79] = 0;
   out_285691373668569915[80] = 1;
}
void h_25(double *state, double *unused, double *out_4174992519817443698) {
   out_4174992519817443698[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2624001866442726016) {
   out_2624001866442726016[0] = 0;
   out_2624001866442726016[1] = 0;
   out_2624001866442726016[2] = 0;
   out_2624001866442726016[3] = 0;
   out_2624001866442726016[4] = 0;
   out_2624001866442726016[5] = 0;
   out_2624001866442726016[6] = 1;
   out_2624001866442726016[7] = 0;
   out_2624001866442726016[8] = 0;
}
void h_24(double *state, double *unused, double *out_9139150395481525125) {
   out_9139150395481525125[0] = state[4];
   out_9139150395481525125[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4801216290049875989) {
   out_4801216290049875989[0] = 0;
   out_4801216290049875989[1] = 0;
   out_4801216290049875989[2] = 0;
   out_4801216290049875989[3] = 0;
   out_4801216290049875989[4] = 1;
   out_4801216290049875989[5] = 0;
   out_4801216290049875989[6] = 0;
   out_4801216290049875989[7] = 0;
   out_4801216290049875989[8] = 0;
   out_4801216290049875989[9] = 0;
   out_4801216290049875989[10] = 0;
   out_4801216290049875989[11] = 0;
   out_4801216290049875989[12] = 0;
   out_4801216290049875989[13] = 0;
   out_4801216290049875989[14] = 1;
   out_4801216290049875989[15] = 0;
   out_4801216290049875989[16] = 0;
   out_4801216290049875989[17] = 0;
}
void h_30(double *state, double *unused, double *out_3588272737003228743) {
   out_3588272737003228743[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5142334824949974643) {
   out_5142334824949974643[0] = 0;
   out_5142334824949974643[1] = 0;
   out_5142334824949974643[2] = 0;
   out_5142334824949974643[3] = 0;
   out_5142334824949974643[4] = 1;
   out_5142334824949974643[5] = 0;
   out_5142334824949974643[6] = 0;
   out_5142334824949974643[7] = 0;
   out_5142334824949974643[8] = 0;
}
void h_26(double *state, double *unused, double *out_551486685558896058) {
   out_551486685558896058[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1117501452431330208) {
   out_1117501452431330208[0] = 0;
   out_1117501452431330208[1] = 0;
   out_1117501452431330208[2] = 0;
   out_1117501452431330208[3] = 0;
   out_1117501452431330208[4] = 0;
   out_1117501452431330208[5] = 0;
   out_1117501452431330208[6] = 0;
   out_1117501452431330208[7] = 1;
   out_1117501452431330208[8] = 0;
}
void h_27(double *state, double *unused, double *out_8169232118158913528) {
   out_8169232118158913528[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2967571513149549732) {
   out_2967571513149549732[0] = 0;
   out_2967571513149549732[1] = 0;
   out_2967571513149549732[2] = 0;
   out_2967571513149549732[3] = 1;
   out_2967571513149549732[4] = 0;
   out_2967571513149549732[5] = 0;
   out_2967571513149549732[6] = 0;
   out_2967571513149549732[7] = 0;
   out_2967571513149549732[8] = 0;
}
void h_29(double *state, double *unused, double *out_1413400530657616902) {
   out_1413400530657616902[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5652566169264366827) {
   out_5652566169264366827[0] = 0;
   out_5652566169264366827[1] = 1;
   out_5652566169264366827[2] = 0;
   out_5652566169264366827[3] = 0;
   out_5652566169264366827[4] = 0;
   out_5652566169264366827[5] = 0;
   out_5652566169264366827[6] = 0;
   out_5652566169264366827[7] = 0;
   out_5652566169264366827[8] = 0;
}
void h_28(double *state, double *unused, double *out_7240725904106305403) {
   out_7240725904106305403[0] = state[0];
}
void H_28(double *state, double *unused, double *out_570167152194836253) {
   out_570167152194836253[0] = 1;
   out_570167152194836253[1] = 0;
   out_570167152194836253[2] = 0;
   out_570167152194836253[3] = 0;
   out_570167152194836253[4] = 0;
   out_570167152194836253[5] = 0;
   out_570167152194836253[6] = 0;
   out_570167152194836253[7] = 0;
   out_570167152194836253[8] = 0;
}
void h_31(double *state, double *unused, double *out_5775090872628367926) {
   out_5775090872628367926[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1743709554664681684) {
   out_1743709554664681684[0] = 0;
   out_1743709554664681684[1] = 0;
   out_1743709554664681684[2] = 0;
   out_1743709554664681684[3] = 0;
   out_1743709554664681684[4] = 0;
   out_1743709554664681684[5] = 0;
   out_1743709554664681684[6] = 0;
   out_1743709554664681684[7] = 0;
   out_1743709554664681684[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6221613173377875958) {
  err_fun(nom_x, delta_x, out_6221613173377875958);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1873114568696790253) {
  inv_err_fun(nom_x, true_x, out_1873114568696790253);
}
void car_H_mod_fun(double *state, double *out_9120391434276283067) {
  H_mod_fun(state, out_9120391434276283067);
}
void car_f_fun(double *state, double dt, double *out_8908813278594628189) {
  f_fun(state,  dt, out_8908813278594628189);
}
void car_F_fun(double *state, double dt, double *out_285691373668569915) {
  F_fun(state,  dt, out_285691373668569915);
}
void car_h_25(double *state, double *unused, double *out_4174992519817443698) {
  h_25(state, unused, out_4174992519817443698);
}
void car_H_25(double *state, double *unused, double *out_2624001866442726016) {
  H_25(state, unused, out_2624001866442726016);
}
void car_h_24(double *state, double *unused, double *out_9139150395481525125) {
  h_24(state, unused, out_9139150395481525125);
}
void car_H_24(double *state, double *unused, double *out_4801216290049875989) {
  H_24(state, unused, out_4801216290049875989);
}
void car_h_30(double *state, double *unused, double *out_3588272737003228743) {
  h_30(state, unused, out_3588272737003228743);
}
void car_H_30(double *state, double *unused, double *out_5142334824949974643) {
  H_30(state, unused, out_5142334824949974643);
}
void car_h_26(double *state, double *unused, double *out_551486685558896058) {
  h_26(state, unused, out_551486685558896058);
}
void car_H_26(double *state, double *unused, double *out_1117501452431330208) {
  H_26(state, unused, out_1117501452431330208);
}
void car_h_27(double *state, double *unused, double *out_8169232118158913528) {
  h_27(state, unused, out_8169232118158913528);
}
void car_H_27(double *state, double *unused, double *out_2967571513149549732) {
  H_27(state, unused, out_2967571513149549732);
}
void car_h_29(double *state, double *unused, double *out_1413400530657616902) {
  h_29(state, unused, out_1413400530657616902);
}
void car_H_29(double *state, double *unused, double *out_5652566169264366827) {
  H_29(state, unused, out_5652566169264366827);
}
void car_h_28(double *state, double *unused, double *out_7240725904106305403) {
  h_28(state, unused, out_7240725904106305403);
}
void car_H_28(double *state, double *unused, double *out_570167152194836253) {
  H_28(state, unused, out_570167152194836253);
}
void car_h_31(double *state, double *unused, double *out_5775090872628367926) {
  h_31(state, unused, out_5775090872628367926);
}
void car_H_31(double *state, double *unused, double *out_1743709554664681684) {
  H_31(state, unused, out_1743709554664681684);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
