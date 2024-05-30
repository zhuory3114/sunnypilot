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
void err_fun(double *nom_x, double *delta_x, double *out_8677092309422431987) {
   out_8677092309422431987[0] = delta_x[0] + nom_x[0];
   out_8677092309422431987[1] = delta_x[1] + nom_x[1];
   out_8677092309422431987[2] = delta_x[2] + nom_x[2];
   out_8677092309422431987[3] = delta_x[3] + nom_x[3];
   out_8677092309422431987[4] = delta_x[4] + nom_x[4];
   out_8677092309422431987[5] = delta_x[5] + nom_x[5];
   out_8677092309422431987[6] = delta_x[6] + nom_x[6];
   out_8677092309422431987[7] = delta_x[7] + nom_x[7];
   out_8677092309422431987[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3536033552698117999) {
   out_3536033552698117999[0] = -nom_x[0] + true_x[0];
   out_3536033552698117999[1] = -nom_x[1] + true_x[1];
   out_3536033552698117999[2] = -nom_x[2] + true_x[2];
   out_3536033552698117999[3] = -nom_x[3] + true_x[3];
   out_3536033552698117999[4] = -nom_x[4] + true_x[4];
   out_3536033552698117999[5] = -nom_x[5] + true_x[5];
   out_3536033552698117999[6] = -nom_x[6] + true_x[6];
   out_3536033552698117999[7] = -nom_x[7] + true_x[7];
   out_3536033552698117999[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8253580055064896978) {
   out_8253580055064896978[0] = 1.0;
   out_8253580055064896978[1] = 0;
   out_8253580055064896978[2] = 0;
   out_8253580055064896978[3] = 0;
   out_8253580055064896978[4] = 0;
   out_8253580055064896978[5] = 0;
   out_8253580055064896978[6] = 0;
   out_8253580055064896978[7] = 0;
   out_8253580055064896978[8] = 0;
   out_8253580055064896978[9] = 0;
   out_8253580055064896978[10] = 1.0;
   out_8253580055064896978[11] = 0;
   out_8253580055064896978[12] = 0;
   out_8253580055064896978[13] = 0;
   out_8253580055064896978[14] = 0;
   out_8253580055064896978[15] = 0;
   out_8253580055064896978[16] = 0;
   out_8253580055064896978[17] = 0;
   out_8253580055064896978[18] = 0;
   out_8253580055064896978[19] = 0;
   out_8253580055064896978[20] = 1.0;
   out_8253580055064896978[21] = 0;
   out_8253580055064896978[22] = 0;
   out_8253580055064896978[23] = 0;
   out_8253580055064896978[24] = 0;
   out_8253580055064896978[25] = 0;
   out_8253580055064896978[26] = 0;
   out_8253580055064896978[27] = 0;
   out_8253580055064896978[28] = 0;
   out_8253580055064896978[29] = 0;
   out_8253580055064896978[30] = 1.0;
   out_8253580055064896978[31] = 0;
   out_8253580055064896978[32] = 0;
   out_8253580055064896978[33] = 0;
   out_8253580055064896978[34] = 0;
   out_8253580055064896978[35] = 0;
   out_8253580055064896978[36] = 0;
   out_8253580055064896978[37] = 0;
   out_8253580055064896978[38] = 0;
   out_8253580055064896978[39] = 0;
   out_8253580055064896978[40] = 1.0;
   out_8253580055064896978[41] = 0;
   out_8253580055064896978[42] = 0;
   out_8253580055064896978[43] = 0;
   out_8253580055064896978[44] = 0;
   out_8253580055064896978[45] = 0;
   out_8253580055064896978[46] = 0;
   out_8253580055064896978[47] = 0;
   out_8253580055064896978[48] = 0;
   out_8253580055064896978[49] = 0;
   out_8253580055064896978[50] = 1.0;
   out_8253580055064896978[51] = 0;
   out_8253580055064896978[52] = 0;
   out_8253580055064896978[53] = 0;
   out_8253580055064896978[54] = 0;
   out_8253580055064896978[55] = 0;
   out_8253580055064896978[56] = 0;
   out_8253580055064896978[57] = 0;
   out_8253580055064896978[58] = 0;
   out_8253580055064896978[59] = 0;
   out_8253580055064896978[60] = 1.0;
   out_8253580055064896978[61] = 0;
   out_8253580055064896978[62] = 0;
   out_8253580055064896978[63] = 0;
   out_8253580055064896978[64] = 0;
   out_8253580055064896978[65] = 0;
   out_8253580055064896978[66] = 0;
   out_8253580055064896978[67] = 0;
   out_8253580055064896978[68] = 0;
   out_8253580055064896978[69] = 0;
   out_8253580055064896978[70] = 1.0;
   out_8253580055064896978[71] = 0;
   out_8253580055064896978[72] = 0;
   out_8253580055064896978[73] = 0;
   out_8253580055064896978[74] = 0;
   out_8253580055064896978[75] = 0;
   out_8253580055064896978[76] = 0;
   out_8253580055064896978[77] = 0;
   out_8253580055064896978[78] = 0;
   out_8253580055064896978[79] = 0;
   out_8253580055064896978[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5986056227021830781) {
   out_5986056227021830781[0] = state[0];
   out_5986056227021830781[1] = state[1];
   out_5986056227021830781[2] = state[2];
   out_5986056227021830781[3] = state[3];
   out_5986056227021830781[4] = state[4];
   out_5986056227021830781[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5986056227021830781[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5986056227021830781[7] = state[7];
   out_5986056227021830781[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8361123820946526770) {
   out_8361123820946526770[0] = 1;
   out_8361123820946526770[1] = 0;
   out_8361123820946526770[2] = 0;
   out_8361123820946526770[3] = 0;
   out_8361123820946526770[4] = 0;
   out_8361123820946526770[5] = 0;
   out_8361123820946526770[6] = 0;
   out_8361123820946526770[7] = 0;
   out_8361123820946526770[8] = 0;
   out_8361123820946526770[9] = 0;
   out_8361123820946526770[10] = 1;
   out_8361123820946526770[11] = 0;
   out_8361123820946526770[12] = 0;
   out_8361123820946526770[13] = 0;
   out_8361123820946526770[14] = 0;
   out_8361123820946526770[15] = 0;
   out_8361123820946526770[16] = 0;
   out_8361123820946526770[17] = 0;
   out_8361123820946526770[18] = 0;
   out_8361123820946526770[19] = 0;
   out_8361123820946526770[20] = 1;
   out_8361123820946526770[21] = 0;
   out_8361123820946526770[22] = 0;
   out_8361123820946526770[23] = 0;
   out_8361123820946526770[24] = 0;
   out_8361123820946526770[25] = 0;
   out_8361123820946526770[26] = 0;
   out_8361123820946526770[27] = 0;
   out_8361123820946526770[28] = 0;
   out_8361123820946526770[29] = 0;
   out_8361123820946526770[30] = 1;
   out_8361123820946526770[31] = 0;
   out_8361123820946526770[32] = 0;
   out_8361123820946526770[33] = 0;
   out_8361123820946526770[34] = 0;
   out_8361123820946526770[35] = 0;
   out_8361123820946526770[36] = 0;
   out_8361123820946526770[37] = 0;
   out_8361123820946526770[38] = 0;
   out_8361123820946526770[39] = 0;
   out_8361123820946526770[40] = 1;
   out_8361123820946526770[41] = 0;
   out_8361123820946526770[42] = 0;
   out_8361123820946526770[43] = 0;
   out_8361123820946526770[44] = 0;
   out_8361123820946526770[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8361123820946526770[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8361123820946526770[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8361123820946526770[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8361123820946526770[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8361123820946526770[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8361123820946526770[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8361123820946526770[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8361123820946526770[53] = -9.8000000000000007*dt;
   out_8361123820946526770[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8361123820946526770[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8361123820946526770[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8361123820946526770[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8361123820946526770[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8361123820946526770[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8361123820946526770[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8361123820946526770[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8361123820946526770[62] = 0;
   out_8361123820946526770[63] = 0;
   out_8361123820946526770[64] = 0;
   out_8361123820946526770[65] = 0;
   out_8361123820946526770[66] = 0;
   out_8361123820946526770[67] = 0;
   out_8361123820946526770[68] = 0;
   out_8361123820946526770[69] = 0;
   out_8361123820946526770[70] = 1;
   out_8361123820946526770[71] = 0;
   out_8361123820946526770[72] = 0;
   out_8361123820946526770[73] = 0;
   out_8361123820946526770[74] = 0;
   out_8361123820946526770[75] = 0;
   out_8361123820946526770[76] = 0;
   out_8361123820946526770[77] = 0;
   out_8361123820946526770[78] = 0;
   out_8361123820946526770[79] = 0;
   out_8361123820946526770[80] = 1;
}
void h_25(double *state, double *unused, double *out_7030946204566009936) {
   out_7030946204566009936[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3735175456469797709) {
   out_3735175456469797709[0] = 0;
   out_3735175456469797709[1] = 0;
   out_3735175456469797709[2] = 0;
   out_3735175456469797709[3] = 0;
   out_3735175456469797709[4] = 0;
   out_3735175456469797709[5] = 0;
   out_3735175456469797709[6] = 1;
   out_3735175456469797709[7] = 0;
   out_3735175456469797709[8] = 0;
}
void h_24(double *state, double *unused, double *out_1397172579223915339) {
   out_1397172579223915339[0] = state[4];
   out_1397172579223915339[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8900182986216976942) {
   out_8900182986216976942[0] = 0;
   out_8900182986216976942[1] = 0;
   out_8900182986216976942[2] = 0;
   out_8900182986216976942[3] = 0;
   out_8900182986216976942[4] = 1;
   out_8900182986216976942[5] = 0;
   out_8900182986216976942[6] = 0;
   out_8900182986216976942[7] = 0;
   out_8900182986216976942[8] = 0;
   out_8900182986216976942[9] = 0;
   out_8900182986216976942[10] = 0;
   out_8900182986216976942[11] = 0;
   out_8900182986216976942[12] = 0;
   out_8900182986216976942[13] = 0;
   out_8900182986216976942[14] = 1;
   out_8900182986216976942[15] = 0;
   out_8900182986216976942[16] = 0;
   out_8900182986216976942[17] = 0;
}
void h_30(double *state, double *unused, double *out_3393496676311987586) {
   out_3393496676311987586[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7794878275748137152) {
   out_7794878275748137152[0] = 0;
   out_7794878275748137152[1] = 0;
   out_7794878275748137152[2] = 0;
   out_7794878275748137152[3] = 0;
   out_7794878275748137152[4] = 1;
   out_7794878275748137152[5] = 0;
   out_7794878275748137152[6] = 0;
   out_7794878275748137152[7] = 0;
   out_7794878275748137152[8] = 0;
}
void h_26(double *state, double *unused, double *out_5588472495694853179) {
   out_5588472495694853179[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7039701426230598310) {
   out_7039701426230598310[0] = 0;
   out_7039701426230598310[1] = 0;
   out_7039701426230598310[2] = 0;
   out_7039701426230598310[3] = 0;
   out_7039701426230598310[4] = 0;
   out_7039701426230598310[5] = 0;
   out_7039701426230598310[6] = 0;
   out_7039701426230598310[7] = 1;
   out_7039701426230598310[8] = 0;
}
void h_27(double *state, double *unused, double *out_2029272936905164291) {
   out_2029272936905164291[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8477102486160989553) {
   out_8477102486160989553[0] = 0;
   out_8477102486160989553[1] = 0;
   out_8477102486160989553[2] = 0;
   out_8477102486160989553[3] = 1;
   out_8477102486160989553[4] = 0;
   out_8477102486160989553[5] = 0;
   out_8477102486160989553[6] = 0;
   out_8477102486160989553[7] = 0;
   out_8477102486160989553[8] = 0;
}
void h_29(double *state, double *unused, double *out_1608176591348858059) {
   out_1608176591348858059[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7284646931433744968) {
   out_7284646931433744968[0] = 0;
   out_7284646931433744968[1] = 1;
   out_7284646931433744968[2] = 0;
   out_7284646931433744968[3] = 0;
   out_7284646931433744968[4] = 0;
   out_7284646931433744968[5] = 0;
   out_7284646931433744968[6] = 0;
   out_7284646931433744968[7] = 0;
   out_7284646931433744968[8] = 0;
}
void h_28(double *state, double *unused, double *out_5669780174728713609) {
   out_5669780174728713609[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6079698125206276074) {
   out_6079698125206276074[0] = 1;
   out_6079698125206276074[1] = 0;
   out_6079698125206276074[2] = 0;
   out_6079698125206276074[3] = 0;
   out_6079698125206276074[4] = 0;
   out_6079698125206276074[5] = 0;
   out_6079698125206276074[6] = 0;
   out_6079698125206276074[7] = 0;
   out_6079698125206276074[8] = 0;
}
void h_31(double *state, double *unused, double *out_7554114392715040028) {
   out_7554114392715040028[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3765821418346758137) {
   out_3765821418346758137[0] = 0;
   out_3765821418346758137[1] = 0;
   out_3765821418346758137[2] = 0;
   out_3765821418346758137[3] = 0;
   out_3765821418346758137[4] = 0;
   out_3765821418346758137[5] = 0;
   out_3765821418346758137[6] = 0;
   out_3765821418346758137[7] = 0;
   out_3765821418346758137[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_8677092309422431987) {
  err_fun(nom_x, delta_x, out_8677092309422431987);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3536033552698117999) {
  inv_err_fun(nom_x, true_x, out_3536033552698117999);
}
void car_H_mod_fun(double *state, double *out_8253580055064896978) {
  H_mod_fun(state, out_8253580055064896978);
}
void car_f_fun(double *state, double dt, double *out_5986056227021830781) {
  f_fun(state,  dt, out_5986056227021830781);
}
void car_F_fun(double *state, double dt, double *out_8361123820946526770) {
  F_fun(state,  dt, out_8361123820946526770);
}
void car_h_25(double *state, double *unused, double *out_7030946204566009936) {
  h_25(state, unused, out_7030946204566009936);
}
void car_H_25(double *state, double *unused, double *out_3735175456469797709) {
  H_25(state, unused, out_3735175456469797709);
}
void car_h_24(double *state, double *unused, double *out_1397172579223915339) {
  h_24(state, unused, out_1397172579223915339);
}
void car_H_24(double *state, double *unused, double *out_8900182986216976942) {
  H_24(state, unused, out_8900182986216976942);
}
void car_h_30(double *state, double *unused, double *out_3393496676311987586) {
  h_30(state, unused, out_3393496676311987586);
}
void car_H_30(double *state, double *unused, double *out_7794878275748137152) {
  H_30(state, unused, out_7794878275748137152);
}
void car_h_26(double *state, double *unused, double *out_5588472495694853179) {
  h_26(state, unused, out_5588472495694853179);
}
void car_H_26(double *state, double *unused, double *out_7039701426230598310) {
  H_26(state, unused, out_7039701426230598310);
}
void car_h_27(double *state, double *unused, double *out_2029272936905164291) {
  h_27(state, unused, out_2029272936905164291);
}
void car_H_27(double *state, double *unused, double *out_8477102486160989553) {
  H_27(state, unused, out_8477102486160989553);
}
void car_h_29(double *state, double *unused, double *out_1608176591348858059) {
  h_29(state, unused, out_1608176591348858059);
}
void car_H_29(double *state, double *unused, double *out_7284646931433744968) {
  H_29(state, unused, out_7284646931433744968);
}
void car_h_28(double *state, double *unused, double *out_5669780174728713609) {
  h_28(state, unused, out_5669780174728713609);
}
void car_H_28(double *state, double *unused, double *out_6079698125206276074) {
  H_28(state, unused, out_6079698125206276074);
}
void car_h_31(double *state, double *unused, double *out_7554114392715040028) {
  h_31(state, unused, out_7554114392715040028);
}
void car_H_31(double *state, double *unused, double *out_3765821418346758137) {
  H_31(state, unused, out_3765821418346758137);
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
