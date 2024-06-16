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
void err_fun(double *nom_x, double *delta_x, double *out_3930116420352682607) {
   out_3930116420352682607[0] = delta_x[0] + nom_x[0];
   out_3930116420352682607[1] = delta_x[1] + nom_x[1];
   out_3930116420352682607[2] = delta_x[2] + nom_x[2];
   out_3930116420352682607[3] = delta_x[3] + nom_x[3];
   out_3930116420352682607[4] = delta_x[4] + nom_x[4];
   out_3930116420352682607[5] = delta_x[5] + nom_x[5];
   out_3930116420352682607[6] = delta_x[6] + nom_x[6];
   out_3930116420352682607[7] = delta_x[7] + nom_x[7];
   out_3930116420352682607[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8332055307244711483) {
   out_8332055307244711483[0] = -nom_x[0] + true_x[0];
   out_8332055307244711483[1] = -nom_x[1] + true_x[1];
   out_8332055307244711483[2] = -nom_x[2] + true_x[2];
   out_8332055307244711483[3] = -nom_x[3] + true_x[3];
   out_8332055307244711483[4] = -nom_x[4] + true_x[4];
   out_8332055307244711483[5] = -nom_x[5] + true_x[5];
   out_8332055307244711483[6] = -nom_x[6] + true_x[6];
   out_8332055307244711483[7] = -nom_x[7] + true_x[7];
   out_8332055307244711483[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1183099394398996416) {
   out_1183099394398996416[0] = 1.0;
   out_1183099394398996416[1] = 0;
   out_1183099394398996416[2] = 0;
   out_1183099394398996416[3] = 0;
   out_1183099394398996416[4] = 0;
   out_1183099394398996416[5] = 0;
   out_1183099394398996416[6] = 0;
   out_1183099394398996416[7] = 0;
   out_1183099394398996416[8] = 0;
   out_1183099394398996416[9] = 0;
   out_1183099394398996416[10] = 1.0;
   out_1183099394398996416[11] = 0;
   out_1183099394398996416[12] = 0;
   out_1183099394398996416[13] = 0;
   out_1183099394398996416[14] = 0;
   out_1183099394398996416[15] = 0;
   out_1183099394398996416[16] = 0;
   out_1183099394398996416[17] = 0;
   out_1183099394398996416[18] = 0;
   out_1183099394398996416[19] = 0;
   out_1183099394398996416[20] = 1.0;
   out_1183099394398996416[21] = 0;
   out_1183099394398996416[22] = 0;
   out_1183099394398996416[23] = 0;
   out_1183099394398996416[24] = 0;
   out_1183099394398996416[25] = 0;
   out_1183099394398996416[26] = 0;
   out_1183099394398996416[27] = 0;
   out_1183099394398996416[28] = 0;
   out_1183099394398996416[29] = 0;
   out_1183099394398996416[30] = 1.0;
   out_1183099394398996416[31] = 0;
   out_1183099394398996416[32] = 0;
   out_1183099394398996416[33] = 0;
   out_1183099394398996416[34] = 0;
   out_1183099394398996416[35] = 0;
   out_1183099394398996416[36] = 0;
   out_1183099394398996416[37] = 0;
   out_1183099394398996416[38] = 0;
   out_1183099394398996416[39] = 0;
   out_1183099394398996416[40] = 1.0;
   out_1183099394398996416[41] = 0;
   out_1183099394398996416[42] = 0;
   out_1183099394398996416[43] = 0;
   out_1183099394398996416[44] = 0;
   out_1183099394398996416[45] = 0;
   out_1183099394398996416[46] = 0;
   out_1183099394398996416[47] = 0;
   out_1183099394398996416[48] = 0;
   out_1183099394398996416[49] = 0;
   out_1183099394398996416[50] = 1.0;
   out_1183099394398996416[51] = 0;
   out_1183099394398996416[52] = 0;
   out_1183099394398996416[53] = 0;
   out_1183099394398996416[54] = 0;
   out_1183099394398996416[55] = 0;
   out_1183099394398996416[56] = 0;
   out_1183099394398996416[57] = 0;
   out_1183099394398996416[58] = 0;
   out_1183099394398996416[59] = 0;
   out_1183099394398996416[60] = 1.0;
   out_1183099394398996416[61] = 0;
   out_1183099394398996416[62] = 0;
   out_1183099394398996416[63] = 0;
   out_1183099394398996416[64] = 0;
   out_1183099394398996416[65] = 0;
   out_1183099394398996416[66] = 0;
   out_1183099394398996416[67] = 0;
   out_1183099394398996416[68] = 0;
   out_1183099394398996416[69] = 0;
   out_1183099394398996416[70] = 1.0;
   out_1183099394398996416[71] = 0;
   out_1183099394398996416[72] = 0;
   out_1183099394398996416[73] = 0;
   out_1183099394398996416[74] = 0;
   out_1183099394398996416[75] = 0;
   out_1183099394398996416[76] = 0;
   out_1183099394398996416[77] = 0;
   out_1183099394398996416[78] = 0;
   out_1183099394398996416[79] = 0;
   out_1183099394398996416[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_3974371446758158566) {
   out_3974371446758158566[0] = state[0];
   out_3974371446758158566[1] = state[1];
   out_3974371446758158566[2] = state[2];
   out_3974371446758158566[3] = state[3];
   out_3974371446758158566[4] = state[4];
   out_3974371446758158566[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3974371446758158566[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3974371446758158566[7] = state[7];
   out_3974371446758158566[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1509814126693874228) {
   out_1509814126693874228[0] = 1;
   out_1509814126693874228[1] = 0;
   out_1509814126693874228[2] = 0;
   out_1509814126693874228[3] = 0;
   out_1509814126693874228[4] = 0;
   out_1509814126693874228[5] = 0;
   out_1509814126693874228[6] = 0;
   out_1509814126693874228[7] = 0;
   out_1509814126693874228[8] = 0;
   out_1509814126693874228[9] = 0;
   out_1509814126693874228[10] = 1;
   out_1509814126693874228[11] = 0;
   out_1509814126693874228[12] = 0;
   out_1509814126693874228[13] = 0;
   out_1509814126693874228[14] = 0;
   out_1509814126693874228[15] = 0;
   out_1509814126693874228[16] = 0;
   out_1509814126693874228[17] = 0;
   out_1509814126693874228[18] = 0;
   out_1509814126693874228[19] = 0;
   out_1509814126693874228[20] = 1;
   out_1509814126693874228[21] = 0;
   out_1509814126693874228[22] = 0;
   out_1509814126693874228[23] = 0;
   out_1509814126693874228[24] = 0;
   out_1509814126693874228[25] = 0;
   out_1509814126693874228[26] = 0;
   out_1509814126693874228[27] = 0;
   out_1509814126693874228[28] = 0;
   out_1509814126693874228[29] = 0;
   out_1509814126693874228[30] = 1;
   out_1509814126693874228[31] = 0;
   out_1509814126693874228[32] = 0;
   out_1509814126693874228[33] = 0;
   out_1509814126693874228[34] = 0;
   out_1509814126693874228[35] = 0;
   out_1509814126693874228[36] = 0;
   out_1509814126693874228[37] = 0;
   out_1509814126693874228[38] = 0;
   out_1509814126693874228[39] = 0;
   out_1509814126693874228[40] = 1;
   out_1509814126693874228[41] = 0;
   out_1509814126693874228[42] = 0;
   out_1509814126693874228[43] = 0;
   out_1509814126693874228[44] = 0;
   out_1509814126693874228[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1509814126693874228[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1509814126693874228[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1509814126693874228[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1509814126693874228[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1509814126693874228[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1509814126693874228[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1509814126693874228[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1509814126693874228[53] = -9.8000000000000007*dt;
   out_1509814126693874228[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1509814126693874228[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1509814126693874228[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1509814126693874228[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1509814126693874228[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1509814126693874228[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1509814126693874228[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1509814126693874228[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1509814126693874228[62] = 0;
   out_1509814126693874228[63] = 0;
   out_1509814126693874228[64] = 0;
   out_1509814126693874228[65] = 0;
   out_1509814126693874228[66] = 0;
   out_1509814126693874228[67] = 0;
   out_1509814126693874228[68] = 0;
   out_1509814126693874228[69] = 0;
   out_1509814126693874228[70] = 1;
   out_1509814126693874228[71] = 0;
   out_1509814126693874228[72] = 0;
   out_1509814126693874228[73] = 0;
   out_1509814126693874228[74] = 0;
   out_1509814126693874228[75] = 0;
   out_1509814126693874228[76] = 0;
   out_1509814126693874228[77] = 0;
   out_1509814126693874228[78] = 0;
   out_1509814126693874228[79] = 0;
   out_1509814126693874228[80] = 1;
}
void h_25(double *state, double *unused, double *out_479281131967255891) {
   out_479281131967255891[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2945433631051225152) {
   out_2945433631051225152[0] = 0;
   out_2945433631051225152[1] = 0;
   out_2945433631051225152[2] = 0;
   out_2945433631051225152[3] = 0;
   out_2945433631051225152[4] = 0;
   out_2945433631051225152[5] = 0;
   out_2945433631051225152[6] = 1;
   out_2945433631051225152[7] = 0;
   out_2945433631051225152[8] = 0;
}
void h_24(double *state, double *unused, double *out_2798532925715571291) {
   out_2798532925715571291[0] = state[4];
   out_2798532925715571291[1] = state[5];
}
void H_24(double *state, double *unused, double *out_772784032045725586) {
   out_772784032045725586[0] = 0;
   out_772784032045725586[1] = 0;
   out_772784032045725586[2] = 0;
   out_772784032045725586[3] = 0;
   out_772784032045725586[4] = 1;
   out_772784032045725586[5] = 0;
   out_772784032045725586[6] = 0;
   out_772784032045725586[7] = 0;
   out_772784032045725586[8] = 0;
   out_772784032045725586[9] = 0;
   out_772784032045725586[10] = 0;
   out_772784032045725586[11] = 0;
   out_772784032045725586[12] = 0;
   out_772784032045725586[13] = 0;
   out_772784032045725586[14] = 1;
   out_772784032045725586[15] = 0;
   out_772784032045725586[16] = 0;
   out_772784032045725586[17] = 0;
}
void h_30(double *state, double *unused, double *out_2517111792866655269) {
   out_2517111792866655269[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5463766589558473779) {
   out_5463766589558473779[0] = 0;
   out_5463766589558473779[1] = 0;
   out_5463766589558473779[2] = 0;
   out_5463766589558473779[3] = 0;
   out_5463766589558473779[4] = 1;
   out_5463766589558473779[5] = 0;
   out_5463766589558473779[6] = 0;
   out_5463766589558473779[7] = 0;
   out_5463766589558473779[8] = 0;
}
void h_26(double *state, double *unused, double *out_6194954447731099144) {
   out_6194954447731099144[0] = state[7];
}
void H_26(double *state, double *unused, double *out_796069687822831072) {
   out_796069687822831072[0] = 0;
   out_796069687822831072[1] = 0;
   out_796069687822831072[2] = 0;
   out_796069687822831072[3] = 0;
   out_796069687822831072[4] = 0;
   out_796069687822831072[5] = 0;
   out_796069687822831072[6] = 0;
   out_796069687822831072[7] = 1;
   out_796069687822831072[8] = 0;
}
void h_27(double *state, double *unused, double *out_1422790984868918326) {
   out_1422790984868918326[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3757026010876807957) {
   out_3757026010876807957[0] = 0;
   out_3757026010876807957[1] = 0;
   out_3757026010876807957[2] = 0;
   out_3757026010876807957[3] = 1;
   out_3757026010876807957[4] = 0;
   out_3757026010876807957[5] = 0;
   out_3757026010876807957[6] = 0;
   out_3757026010876807957[7] = 0;
   out_3757026010876807957[8] = 0;
}
void h_29(double *state, double *unused, double *out_8159841663947612104) {
   out_8159841663947612104[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1072031354761990862) {
   out_1072031354761990862[0] = 0;
   out_1072031354761990862[1] = 1;
   out_1072031354761990862[2] = 0;
   out_1072031354761990862[3] = 0;
   out_1072031354761990862[4] = 0;
   out_1072031354761990862[5] = 0;
   out_1072031354761990862[6] = 0;
   out_1072031354761990862[7] = 0;
   out_1072031354761990862[8] = 0;
}
void h_28(double *state, double *unused, double *out_5704354826112578834) {
   out_5704354826112578834[0] = state[0];
}
void H_28(double *state, double *unused, double *out_891598916803335389) {
   out_891598916803335389[0] = 1;
   out_891598916803335389[1] = 0;
   out_891598916803335389[2] = 0;
   out_891598916803335389[3] = 0;
   out_891598916803335389[4] = 0;
   out_891598916803335389[5] = 0;
   out_891598916803335389[6] = 0;
   out_891598916803335389[7] = 0;
   out_891598916803335389[8] = 0;
}
void h_31(double *state, double *unused, double *out_8360553782721305076) {
   out_8360553782721305076[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1422277790056182548) {
   out_1422277790056182548[0] = 0;
   out_1422277790056182548[1] = 0;
   out_1422277790056182548[2] = 0;
   out_1422277790056182548[3] = 0;
   out_1422277790056182548[4] = 0;
   out_1422277790056182548[5] = 0;
   out_1422277790056182548[6] = 0;
   out_1422277790056182548[7] = 0;
   out_1422277790056182548[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_3930116420352682607) {
  err_fun(nom_x, delta_x, out_3930116420352682607);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8332055307244711483) {
  inv_err_fun(nom_x, true_x, out_8332055307244711483);
}
void car_H_mod_fun(double *state, double *out_1183099394398996416) {
  H_mod_fun(state, out_1183099394398996416);
}
void car_f_fun(double *state, double dt, double *out_3974371446758158566) {
  f_fun(state,  dt, out_3974371446758158566);
}
void car_F_fun(double *state, double dt, double *out_1509814126693874228) {
  F_fun(state,  dt, out_1509814126693874228);
}
void car_h_25(double *state, double *unused, double *out_479281131967255891) {
  h_25(state, unused, out_479281131967255891);
}
void car_H_25(double *state, double *unused, double *out_2945433631051225152) {
  H_25(state, unused, out_2945433631051225152);
}
void car_h_24(double *state, double *unused, double *out_2798532925715571291) {
  h_24(state, unused, out_2798532925715571291);
}
void car_H_24(double *state, double *unused, double *out_772784032045725586) {
  H_24(state, unused, out_772784032045725586);
}
void car_h_30(double *state, double *unused, double *out_2517111792866655269) {
  h_30(state, unused, out_2517111792866655269);
}
void car_H_30(double *state, double *unused, double *out_5463766589558473779) {
  H_30(state, unused, out_5463766589558473779);
}
void car_h_26(double *state, double *unused, double *out_6194954447731099144) {
  h_26(state, unused, out_6194954447731099144);
}
void car_H_26(double *state, double *unused, double *out_796069687822831072) {
  H_26(state, unused, out_796069687822831072);
}
void car_h_27(double *state, double *unused, double *out_1422790984868918326) {
  h_27(state, unused, out_1422790984868918326);
}
void car_H_27(double *state, double *unused, double *out_3757026010876807957) {
  H_27(state, unused, out_3757026010876807957);
}
void car_h_29(double *state, double *unused, double *out_8159841663947612104) {
  h_29(state, unused, out_8159841663947612104);
}
void car_H_29(double *state, double *unused, double *out_1072031354761990862) {
  H_29(state, unused, out_1072031354761990862);
}
void car_h_28(double *state, double *unused, double *out_5704354826112578834) {
  h_28(state, unused, out_5704354826112578834);
}
void car_H_28(double *state, double *unused, double *out_891598916803335389) {
  H_28(state, unused, out_891598916803335389);
}
void car_h_31(double *state, double *unused, double *out_8360553782721305076) {
  h_31(state, unused, out_8360553782721305076);
}
void car_H_31(double *state, double *unused, double *out_1422277790056182548) {
  H_31(state, unused, out_1422277790056182548);
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
