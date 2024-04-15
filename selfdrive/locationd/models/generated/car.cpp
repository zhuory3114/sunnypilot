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
void err_fun(double *nom_x, double *delta_x, double *out_2907344304532562836) {
   out_2907344304532562836[0] = delta_x[0] + nom_x[0];
   out_2907344304532562836[1] = delta_x[1] + nom_x[1];
   out_2907344304532562836[2] = delta_x[2] + nom_x[2];
   out_2907344304532562836[3] = delta_x[3] + nom_x[3];
   out_2907344304532562836[4] = delta_x[4] + nom_x[4];
   out_2907344304532562836[5] = delta_x[5] + nom_x[5];
   out_2907344304532562836[6] = delta_x[6] + nom_x[6];
   out_2907344304532562836[7] = delta_x[7] + nom_x[7];
   out_2907344304532562836[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8954218720378571779) {
   out_8954218720378571779[0] = -nom_x[0] + true_x[0];
   out_8954218720378571779[1] = -nom_x[1] + true_x[1];
   out_8954218720378571779[2] = -nom_x[2] + true_x[2];
   out_8954218720378571779[3] = -nom_x[3] + true_x[3];
   out_8954218720378571779[4] = -nom_x[4] + true_x[4];
   out_8954218720378571779[5] = -nom_x[5] + true_x[5];
   out_8954218720378571779[6] = -nom_x[6] + true_x[6];
   out_8954218720378571779[7] = -nom_x[7] + true_x[7];
   out_8954218720378571779[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6798369121237450864) {
   out_6798369121237450864[0] = 1.0;
   out_6798369121237450864[1] = 0;
   out_6798369121237450864[2] = 0;
   out_6798369121237450864[3] = 0;
   out_6798369121237450864[4] = 0;
   out_6798369121237450864[5] = 0;
   out_6798369121237450864[6] = 0;
   out_6798369121237450864[7] = 0;
   out_6798369121237450864[8] = 0;
   out_6798369121237450864[9] = 0;
   out_6798369121237450864[10] = 1.0;
   out_6798369121237450864[11] = 0;
   out_6798369121237450864[12] = 0;
   out_6798369121237450864[13] = 0;
   out_6798369121237450864[14] = 0;
   out_6798369121237450864[15] = 0;
   out_6798369121237450864[16] = 0;
   out_6798369121237450864[17] = 0;
   out_6798369121237450864[18] = 0;
   out_6798369121237450864[19] = 0;
   out_6798369121237450864[20] = 1.0;
   out_6798369121237450864[21] = 0;
   out_6798369121237450864[22] = 0;
   out_6798369121237450864[23] = 0;
   out_6798369121237450864[24] = 0;
   out_6798369121237450864[25] = 0;
   out_6798369121237450864[26] = 0;
   out_6798369121237450864[27] = 0;
   out_6798369121237450864[28] = 0;
   out_6798369121237450864[29] = 0;
   out_6798369121237450864[30] = 1.0;
   out_6798369121237450864[31] = 0;
   out_6798369121237450864[32] = 0;
   out_6798369121237450864[33] = 0;
   out_6798369121237450864[34] = 0;
   out_6798369121237450864[35] = 0;
   out_6798369121237450864[36] = 0;
   out_6798369121237450864[37] = 0;
   out_6798369121237450864[38] = 0;
   out_6798369121237450864[39] = 0;
   out_6798369121237450864[40] = 1.0;
   out_6798369121237450864[41] = 0;
   out_6798369121237450864[42] = 0;
   out_6798369121237450864[43] = 0;
   out_6798369121237450864[44] = 0;
   out_6798369121237450864[45] = 0;
   out_6798369121237450864[46] = 0;
   out_6798369121237450864[47] = 0;
   out_6798369121237450864[48] = 0;
   out_6798369121237450864[49] = 0;
   out_6798369121237450864[50] = 1.0;
   out_6798369121237450864[51] = 0;
   out_6798369121237450864[52] = 0;
   out_6798369121237450864[53] = 0;
   out_6798369121237450864[54] = 0;
   out_6798369121237450864[55] = 0;
   out_6798369121237450864[56] = 0;
   out_6798369121237450864[57] = 0;
   out_6798369121237450864[58] = 0;
   out_6798369121237450864[59] = 0;
   out_6798369121237450864[60] = 1.0;
   out_6798369121237450864[61] = 0;
   out_6798369121237450864[62] = 0;
   out_6798369121237450864[63] = 0;
   out_6798369121237450864[64] = 0;
   out_6798369121237450864[65] = 0;
   out_6798369121237450864[66] = 0;
   out_6798369121237450864[67] = 0;
   out_6798369121237450864[68] = 0;
   out_6798369121237450864[69] = 0;
   out_6798369121237450864[70] = 1.0;
   out_6798369121237450864[71] = 0;
   out_6798369121237450864[72] = 0;
   out_6798369121237450864[73] = 0;
   out_6798369121237450864[74] = 0;
   out_6798369121237450864[75] = 0;
   out_6798369121237450864[76] = 0;
   out_6798369121237450864[77] = 0;
   out_6798369121237450864[78] = 0;
   out_6798369121237450864[79] = 0;
   out_6798369121237450864[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6734862407766948837) {
   out_6734862407766948837[0] = state[0];
   out_6734862407766948837[1] = state[1];
   out_6734862407766948837[2] = state[2];
   out_6734862407766948837[3] = state[3];
   out_6734862407766948837[4] = state[4];
   out_6734862407766948837[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6734862407766948837[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6734862407766948837[7] = state[7];
   out_6734862407766948837[8] = state[8];
}
void F_fun(double *state, double dt, double *out_7156118317797004448) {
   out_7156118317797004448[0] = 1;
   out_7156118317797004448[1] = 0;
   out_7156118317797004448[2] = 0;
   out_7156118317797004448[3] = 0;
   out_7156118317797004448[4] = 0;
   out_7156118317797004448[5] = 0;
   out_7156118317797004448[6] = 0;
   out_7156118317797004448[7] = 0;
   out_7156118317797004448[8] = 0;
   out_7156118317797004448[9] = 0;
   out_7156118317797004448[10] = 1;
   out_7156118317797004448[11] = 0;
   out_7156118317797004448[12] = 0;
   out_7156118317797004448[13] = 0;
   out_7156118317797004448[14] = 0;
   out_7156118317797004448[15] = 0;
   out_7156118317797004448[16] = 0;
   out_7156118317797004448[17] = 0;
   out_7156118317797004448[18] = 0;
   out_7156118317797004448[19] = 0;
   out_7156118317797004448[20] = 1;
   out_7156118317797004448[21] = 0;
   out_7156118317797004448[22] = 0;
   out_7156118317797004448[23] = 0;
   out_7156118317797004448[24] = 0;
   out_7156118317797004448[25] = 0;
   out_7156118317797004448[26] = 0;
   out_7156118317797004448[27] = 0;
   out_7156118317797004448[28] = 0;
   out_7156118317797004448[29] = 0;
   out_7156118317797004448[30] = 1;
   out_7156118317797004448[31] = 0;
   out_7156118317797004448[32] = 0;
   out_7156118317797004448[33] = 0;
   out_7156118317797004448[34] = 0;
   out_7156118317797004448[35] = 0;
   out_7156118317797004448[36] = 0;
   out_7156118317797004448[37] = 0;
   out_7156118317797004448[38] = 0;
   out_7156118317797004448[39] = 0;
   out_7156118317797004448[40] = 1;
   out_7156118317797004448[41] = 0;
   out_7156118317797004448[42] = 0;
   out_7156118317797004448[43] = 0;
   out_7156118317797004448[44] = 0;
   out_7156118317797004448[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7156118317797004448[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7156118317797004448[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7156118317797004448[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7156118317797004448[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7156118317797004448[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7156118317797004448[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7156118317797004448[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7156118317797004448[53] = -9.8000000000000007*dt;
   out_7156118317797004448[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7156118317797004448[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7156118317797004448[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7156118317797004448[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7156118317797004448[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7156118317797004448[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7156118317797004448[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7156118317797004448[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7156118317797004448[62] = 0;
   out_7156118317797004448[63] = 0;
   out_7156118317797004448[64] = 0;
   out_7156118317797004448[65] = 0;
   out_7156118317797004448[66] = 0;
   out_7156118317797004448[67] = 0;
   out_7156118317797004448[68] = 0;
   out_7156118317797004448[69] = 0;
   out_7156118317797004448[70] = 1;
   out_7156118317797004448[71] = 0;
   out_7156118317797004448[72] = 0;
   out_7156118317797004448[73] = 0;
   out_7156118317797004448[74] = 0;
   out_7156118317797004448[75] = 0;
   out_7156118317797004448[76] = 0;
   out_7156118317797004448[77] = 0;
   out_7156118317797004448[78] = 0;
   out_7156118317797004448[79] = 0;
   out_7156118317797004448[80] = 1;
}
void h_25(double *state, double *unused, double *out_5350053244966229021) {
   out_5350053244966229021[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2711664023153060623) {
   out_2711664023153060623[0] = 0;
   out_2711664023153060623[1] = 0;
   out_2711664023153060623[2] = 0;
   out_2711664023153060623[3] = 0;
   out_2711664023153060623[4] = 0;
   out_2711664023153060623[5] = 0;
   out_2711664023153060623[6] = 1;
   out_2711664023153060623[7] = 0;
   out_2711664023153060623[8] = 0;
}
void h_24(double *state, double *unused, double *out_5232067606711025854) {
   out_5232067606711025854[0] = state[4];
   out_5232067606711025854[1] = state[5];
}
void H_24(double *state, double *unused, double *out_539014424147561057) {
   out_539014424147561057[0] = 0;
   out_539014424147561057[1] = 0;
   out_539014424147561057[2] = 0;
   out_539014424147561057[3] = 0;
   out_539014424147561057[4] = 1;
   out_539014424147561057[5] = 0;
   out_539014424147561057[6] = 0;
   out_539014424147561057[7] = 0;
   out_539014424147561057[8] = 0;
   out_539014424147561057[9] = 0;
   out_539014424147561057[10] = 0;
   out_539014424147561057[11] = 0;
   out_539014424147561057[12] = 0;
   out_539014424147561057[13] = 0;
   out_539014424147561057[14] = 1;
   out_539014424147561057[15] = 0;
   out_539014424147561057[16] = 0;
   out_539014424147561057[17] = 0;
}
void h_30(double *state, double *unused, double *out_5636132064494434270) {
   out_5636132064494434270[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5229996981660309250) {
   out_5229996981660309250[0] = 0;
   out_5229996981660309250[1] = 0;
   out_5229996981660309250[2] = 0;
   out_5229996981660309250[3] = 0;
   out_5229996981660309250[4] = 1;
   out_5229996981660309250[5] = 0;
   out_5229996981660309250[6] = 0;
   out_5229996981660309250[7] = 0;
   out_5229996981660309250[8] = 0;
}
void h_26(double *state, double *unused, double *out_6554787287256070092) {
   out_6554787287256070092[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1029839295720995601) {
   out_1029839295720995601[0] = 0;
   out_1029839295720995601[1] = 0;
   out_1029839295720995601[2] = 0;
   out_1029839295720995601[3] = 0;
   out_1029839295720995601[4] = 0;
   out_1029839295720995601[5] = 0;
   out_1029839295720995601[6] = 0;
   out_1029839295720995601[7] = 1;
   out_1029839295720995601[8] = 0;
}
void h_27(double *state, double *unused, double *out_3670490167474247421) {
   out_3670490167474247421[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3055233669859884339) {
   out_3055233669859884339[0] = 0;
   out_3055233669859884339[1] = 0;
   out_3055233669859884339[2] = 0;
   out_3055233669859884339[3] = 1;
   out_3055233669859884339[4] = 0;
   out_3055233669859884339[5] = 0;
   out_3055233669859884339[6] = 0;
   out_3055233669859884339[7] = 0;
   out_3055233669859884339[8] = 0;
}
void h_29(double *state, double *unused, double *out_3513303499123118276) {
   out_3513303499123118276[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5740228325974701434) {
   out_5740228325974701434[0] = 0;
   out_5740228325974701434[1] = 1;
   out_5740228325974701434[2] = 0;
   out_5740228325974701434[3] = 0;
   out_5740228325974701434[4] = 0;
   out_5740228325974701434[5] = 0;
   out_5740228325974701434[6] = 0;
   out_5740228325974701434[7] = 0;
   out_5740228325974701434[8] = 0;
}
void h_28(double *state, double *unused, double *out_2728551349939601725) {
   out_2728551349939601725[0] = state[0];
}
void H_28(double *state, double *unused, double *out_657829308905170860) {
   out_657829308905170860[0] = 1;
   out_657829308905170860[1] = 0;
   out_657829308905170860[2] = 0;
   out_657829308905170860[3] = 0;
   out_657829308905170860[4] = 0;
   out_657829308905170860[5] = 0;
   out_657829308905170860[6] = 0;
   out_657829308905170860[7] = 0;
   out_657829308905170860[8] = 0;
}
void h_31(double *state, double *unused, double *out_8987502773220251371) {
   out_8987502773220251371[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1656047397954347077) {
   out_1656047397954347077[0] = 0;
   out_1656047397954347077[1] = 0;
   out_1656047397954347077[2] = 0;
   out_1656047397954347077[3] = 0;
   out_1656047397954347077[4] = 0;
   out_1656047397954347077[5] = 0;
   out_1656047397954347077[6] = 0;
   out_1656047397954347077[7] = 0;
   out_1656047397954347077[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_2907344304532562836) {
  err_fun(nom_x, delta_x, out_2907344304532562836);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8954218720378571779) {
  inv_err_fun(nom_x, true_x, out_8954218720378571779);
}
void car_H_mod_fun(double *state, double *out_6798369121237450864) {
  H_mod_fun(state, out_6798369121237450864);
}
void car_f_fun(double *state, double dt, double *out_6734862407766948837) {
  f_fun(state,  dt, out_6734862407766948837);
}
void car_F_fun(double *state, double dt, double *out_7156118317797004448) {
  F_fun(state,  dt, out_7156118317797004448);
}
void car_h_25(double *state, double *unused, double *out_5350053244966229021) {
  h_25(state, unused, out_5350053244966229021);
}
void car_H_25(double *state, double *unused, double *out_2711664023153060623) {
  H_25(state, unused, out_2711664023153060623);
}
void car_h_24(double *state, double *unused, double *out_5232067606711025854) {
  h_24(state, unused, out_5232067606711025854);
}
void car_H_24(double *state, double *unused, double *out_539014424147561057) {
  H_24(state, unused, out_539014424147561057);
}
void car_h_30(double *state, double *unused, double *out_5636132064494434270) {
  h_30(state, unused, out_5636132064494434270);
}
void car_H_30(double *state, double *unused, double *out_5229996981660309250) {
  H_30(state, unused, out_5229996981660309250);
}
void car_h_26(double *state, double *unused, double *out_6554787287256070092) {
  h_26(state, unused, out_6554787287256070092);
}
void car_H_26(double *state, double *unused, double *out_1029839295720995601) {
  H_26(state, unused, out_1029839295720995601);
}
void car_h_27(double *state, double *unused, double *out_3670490167474247421) {
  h_27(state, unused, out_3670490167474247421);
}
void car_H_27(double *state, double *unused, double *out_3055233669859884339) {
  H_27(state, unused, out_3055233669859884339);
}
void car_h_29(double *state, double *unused, double *out_3513303499123118276) {
  h_29(state, unused, out_3513303499123118276);
}
void car_H_29(double *state, double *unused, double *out_5740228325974701434) {
  H_29(state, unused, out_5740228325974701434);
}
void car_h_28(double *state, double *unused, double *out_2728551349939601725) {
  h_28(state, unused, out_2728551349939601725);
}
void car_H_28(double *state, double *unused, double *out_657829308905170860) {
  H_28(state, unused, out_657829308905170860);
}
void car_h_31(double *state, double *unused, double *out_8987502773220251371) {
  h_31(state, unused, out_8987502773220251371);
}
void car_H_31(double *state, double *unused, double *out_1656047397954347077) {
  H_31(state, unused, out_1656047397954347077);
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
