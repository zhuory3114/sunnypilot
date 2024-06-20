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
void err_fun(double *nom_x, double *delta_x, double *out_6387681849755696540) {
   out_6387681849755696540[0] = delta_x[0] + nom_x[0];
   out_6387681849755696540[1] = delta_x[1] + nom_x[1];
   out_6387681849755696540[2] = delta_x[2] + nom_x[2];
   out_6387681849755696540[3] = delta_x[3] + nom_x[3];
   out_6387681849755696540[4] = delta_x[4] + nom_x[4];
   out_6387681849755696540[5] = delta_x[5] + nom_x[5];
   out_6387681849755696540[6] = delta_x[6] + nom_x[6];
   out_6387681849755696540[7] = delta_x[7] + nom_x[7];
   out_6387681849755696540[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1447910342006919488) {
   out_1447910342006919488[0] = -nom_x[0] + true_x[0];
   out_1447910342006919488[1] = -nom_x[1] + true_x[1];
   out_1447910342006919488[2] = -nom_x[2] + true_x[2];
   out_1447910342006919488[3] = -nom_x[3] + true_x[3];
   out_1447910342006919488[4] = -nom_x[4] + true_x[4];
   out_1447910342006919488[5] = -nom_x[5] + true_x[5];
   out_1447910342006919488[6] = -nom_x[6] + true_x[6];
   out_1447910342006919488[7] = -nom_x[7] + true_x[7];
   out_1447910342006919488[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3269785016151077926) {
   out_3269785016151077926[0] = 1.0;
   out_3269785016151077926[1] = 0;
   out_3269785016151077926[2] = 0;
   out_3269785016151077926[3] = 0;
   out_3269785016151077926[4] = 0;
   out_3269785016151077926[5] = 0;
   out_3269785016151077926[6] = 0;
   out_3269785016151077926[7] = 0;
   out_3269785016151077926[8] = 0;
   out_3269785016151077926[9] = 0;
   out_3269785016151077926[10] = 1.0;
   out_3269785016151077926[11] = 0;
   out_3269785016151077926[12] = 0;
   out_3269785016151077926[13] = 0;
   out_3269785016151077926[14] = 0;
   out_3269785016151077926[15] = 0;
   out_3269785016151077926[16] = 0;
   out_3269785016151077926[17] = 0;
   out_3269785016151077926[18] = 0;
   out_3269785016151077926[19] = 0;
   out_3269785016151077926[20] = 1.0;
   out_3269785016151077926[21] = 0;
   out_3269785016151077926[22] = 0;
   out_3269785016151077926[23] = 0;
   out_3269785016151077926[24] = 0;
   out_3269785016151077926[25] = 0;
   out_3269785016151077926[26] = 0;
   out_3269785016151077926[27] = 0;
   out_3269785016151077926[28] = 0;
   out_3269785016151077926[29] = 0;
   out_3269785016151077926[30] = 1.0;
   out_3269785016151077926[31] = 0;
   out_3269785016151077926[32] = 0;
   out_3269785016151077926[33] = 0;
   out_3269785016151077926[34] = 0;
   out_3269785016151077926[35] = 0;
   out_3269785016151077926[36] = 0;
   out_3269785016151077926[37] = 0;
   out_3269785016151077926[38] = 0;
   out_3269785016151077926[39] = 0;
   out_3269785016151077926[40] = 1.0;
   out_3269785016151077926[41] = 0;
   out_3269785016151077926[42] = 0;
   out_3269785016151077926[43] = 0;
   out_3269785016151077926[44] = 0;
   out_3269785016151077926[45] = 0;
   out_3269785016151077926[46] = 0;
   out_3269785016151077926[47] = 0;
   out_3269785016151077926[48] = 0;
   out_3269785016151077926[49] = 0;
   out_3269785016151077926[50] = 1.0;
   out_3269785016151077926[51] = 0;
   out_3269785016151077926[52] = 0;
   out_3269785016151077926[53] = 0;
   out_3269785016151077926[54] = 0;
   out_3269785016151077926[55] = 0;
   out_3269785016151077926[56] = 0;
   out_3269785016151077926[57] = 0;
   out_3269785016151077926[58] = 0;
   out_3269785016151077926[59] = 0;
   out_3269785016151077926[60] = 1.0;
   out_3269785016151077926[61] = 0;
   out_3269785016151077926[62] = 0;
   out_3269785016151077926[63] = 0;
   out_3269785016151077926[64] = 0;
   out_3269785016151077926[65] = 0;
   out_3269785016151077926[66] = 0;
   out_3269785016151077926[67] = 0;
   out_3269785016151077926[68] = 0;
   out_3269785016151077926[69] = 0;
   out_3269785016151077926[70] = 1.0;
   out_3269785016151077926[71] = 0;
   out_3269785016151077926[72] = 0;
   out_3269785016151077926[73] = 0;
   out_3269785016151077926[74] = 0;
   out_3269785016151077926[75] = 0;
   out_3269785016151077926[76] = 0;
   out_3269785016151077926[77] = 0;
   out_3269785016151077926[78] = 0;
   out_3269785016151077926[79] = 0;
   out_3269785016151077926[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6160736893164847577) {
   out_6160736893164847577[0] = state[0];
   out_6160736893164847577[1] = state[1];
   out_6160736893164847577[2] = state[2];
   out_6160736893164847577[3] = state[3];
   out_6160736893164847577[4] = state[4];
   out_6160736893164847577[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6160736893164847577[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6160736893164847577[7] = state[7];
   out_6160736893164847577[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8208644078465651097) {
   out_8208644078465651097[0] = 1;
   out_8208644078465651097[1] = 0;
   out_8208644078465651097[2] = 0;
   out_8208644078465651097[3] = 0;
   out_8208644078465651097[4] = 0;
   out_8208644078465651097[5] = 0;
   out_8208644078465651097[6] = 0;
   out_8208644078465651097[7] = 0;
   out_8208644078465651097[8] = 0;
   out_8208644078465651097[9] = 0;
   out_8208644078465651097[10] = 1;
   out_8208644078465651097[11] = 0;
   out_8208644078465651097[12] = 0;
   out_8208644078465651097[13] = 0;
   out_8208644078465651097[14] = 0;
   out_8208644078465651097[15] = 0;
   out_8208644078465651097[16] = 0;
   out_8208644078465651097[17] = 0;
   out_8208644078465651097[18] = 0;
   out_8208644078465651097[19] = 0;
   out_8208644078465651097[20] = 1;
   out_8208644078465651097[21] = 0;
   out_8208644078465651097[22] = 0;
   out_8208644078465651097[23] = 0;
   out_8208644078465651097[24] = 0;
   out_8208644078465651097[25] = 0;
   out_8208644078465651097[26] = 0;
   out_8208644078465651097[27] = 0;
   out_8208644078465651097[28] = 0;
   out_8208644078465651097[29] = 0;
   out_8208644078465651097[30] = 1;
   out_8208644078465651097[31] = 0;
   out_8208644078465651097[32] = 0;
   out_8208644078465651097[33] = 0;
   out_8208644078465651097[34] = 0;
   out_8208644078465651097[35] = 0;
   out_8208644078465651097[36] = 0;
   out_8208644078465651097[37] = 0;
   out_8208644078465651097[38] = 0;
   out_8208644078465651097[39] = 0;
   out_8208644078465651097[40] = 1;
   out_8208644078465651097[41] = 0;
   out_8208644078465651097[42] = 0;
   out_8208644078465651097[43] = 0;
   out_8208644078465651097[44] = 0;
   out_8208644078465651097[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8208644078465651097[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8208644078465651097[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8208644078465651097[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8208644078465651097[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8208644078465651097[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8208644078465651097[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8208644078465651097[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8208644078465651097[53] = -9.8000000000000007*dt;
   out_8208644078465651097[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8208644078465651097[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8208644078465651097[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8208644078465651097[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8208644078465651097[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8208644078465651097[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8208644078465651097[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8208644078465651097[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8208644078465651097[62] = 0;
   out_8208644078465651097[63] = 0;
   out_8208644078465651097[64] = 0;
   out_8208644078465651097[65] = 0;
   out_8208644078465651097[66] = 0;
   out_8208644078465651097[67] = 0;
   out_8208644078465651097[68] = 0;
   out_8208644078465651097[69] = 0;
   out_8208644078465651097[70] = 1;
   out_8208644078465651097[71] = 0;
   out_8208644078465651097[72] = 0;
   out_8208644078465651097[73] = 0;
   out_8208644078465651097[74] = 0;
   out_8208644078465651097[75] = 0;
   out_8208644078465651097[76] = 0;
   out_8208644078465651097[77] = 0;
   out_8208644078465651097[78] = 0;
   out_8208644078465651097[79] = 0;
   out_8208644078465651097[80] = 1;
}
void h_25(double *state, double *unused, double *out_4728890029886282683) {
   out_4728890029886282683[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3148005063319919381) {
   out_3148005063319919381[0] = 0;
   out_3148005063319919381[1] = 0;
   out_3148005063319919381[2] = 0;
   out_3148005063319919381[3] = 0;
   out_3148005063319919381[4] = 0;
   out_3148005063319919381[5] = 0;
   out_3148005063319919381[6] = 1;
   out_3148005063319919381[7] = 0;
   out_3148005063319919381[8] = 0;
}
void h_24(double *state, double *unused, double *out_6821372846868724365) {
   out_6821372846868724365[0] = state[4];
   out_6821372846868724365[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3993001672103614780) {
   out_3993001672103614780[0] = 0;
   out_3993001672103614780[1] = 0;
   out_3993001672103614780[2] = 0;
   out_3993001672103614780[3] = 0;
   out_3993001672103614780[4] = 1;
   out_3993001672103614780[5] = 0;
   out_3993001672103614780[6] = 0;
   out_3993001672103614780[7] = 0;
   out_3993001672103614780[8] = 0;
   out_3993001672103614780[9] = 0;
   out_3993001672103614780[10] = 0;
   out_3993001672103614780[11] = 0;
   out_3993001672103614780[12] = 0;
   out_3993001672103614780[13] = 0;
   out_3993001672103614780[14] = 1;
   out_3993001672103614780[15] = 0;
   out_3993001672103614780[16] = 0;
   out_3993001672103614780[17] = 0;
}
void h_30(double *state, double *unused, double *out_2090612130942494452) {
   out_2090612130942494452[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8382048668898015480) {
   out_8382048668898015480[0] = 0;
   out_8382048668898015480[1] = 0;
   out_8382048668898015480[2] = 0;
   out_8382048668898015480[3] = 0;
   out_8382048668898015480[4] = 1;
   out_8382048668898015480[5] = 0;
   out_8382048668898015480[6] = 0;
   out_8382048668898015480[7] = 0;
   out_8382048668898015480[8] = 0;
}
void h_26(double *state, double *unused, double *out_233144791967815326) {
   out_233144791967815326[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6452531033080719982) {
   out_6452531033080719982[0] = 0;
   out_6452531033080719982[1] = 0;
   out_6452531033080719982[2] = 0;
   out_6452531033080719982[3] = 0;
   out_6452531033080719982[4] = 0;
   out_6452531033080719982[5] = 0;
   out_6452531033080719982[6] = 0;
   out_6452531033080719982[7] = 1;
   out_6452531033080719982[8] = 0;
}
void h_27(double *state, double *unused, double *out_7850890224567832796) {
   out_7850890224567832796[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7889932093011111225) {
   out_7889932093011111225[0] = 0;
   out_7889932093011111225[1] = 0;
   out_7889932093011111225[2] = 0;
   out_7889932093011111225[3] = 1;
   out_7889932093011111225[4] = 0;
   out_7889932093011111225[5] = 0;
   out_7889932093011111225[6] = 0;
   out_7889932093011111225[7] = 0;
   out_7889932093011111225[8] = 0;
}
void h_29(double *state, double *unused, double *out_1587810802521718556) {
   out_1587810802521718556[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7871817324583623296) {
   out_7871817324583623296[0] = 0;
   out_7871817324583623296[1] = 1;
   out_7871817324583623296[2] = 0;
   out_7871817324583623296[3] = 0;
   out_7871817324583623296[4] = 0;
   out_7871817324583623296[5] = 0;
   out_7871817324583623296[6] = 0;
   out_7871817324583623296[7] = 0;
   out_7871817324583623296[8] = 0;
}
void h_28(double *state, double *unused, double *out_7991448528333021169) {
   out_7991448528333021169[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1094170349072029618) {
   out_1094170349072029618[0] = 1;
   out_1094170349072029618[1] = 0;
   out_1094170349072029618[2] = 0;
   out_1094170349072029618[3] = 0;
   out_1094170349072029618[4] = 0;
   out_1094170349072029618[5] = 0;
   out_1094170349072029618[6] = 0;
   out_1094170349072029618[7] = 0;
   out_1094170349072029618[8] = 0;
}
void h_31(double *state, double *unused, double *out_1732497105052371523) {
   out_1732497105052371523[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3178651025196879809) {
   out_3178651025196879809[0] = 0;
   out_3178651025196879809[1] = 0;
   out_3178651025196879809[2] = 0;
   out_3178651025196879809[3] = 0;
   out_3178651025196879809[4] = 0;
   out_3178651025196879809[5] = 0;
   out_3178651025196879809[6] = 0;
   out_3178651025196879809[7] = 0;
   out_3178651025196879809[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_6387681849755696540) {
  err_fun(nom_x, delta_x, out_6387681849755696540);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1447910342006919488) {
  inv_err_fun(nom_x, true_x, out_1447910342006919488);
}
void car_H_mod_fun(double *state, double *out_3269785016151077926) {
  H_mod_fun(state, out_3269785016151077926);
}
void car_f_fun(double *state, double dt, double *out_6160736893164847577) {
  f_fun(state,  dt, out_6160736893164847577);
}
void car_F_fun(double *state, double dt, double *out_8208644078465651097) {
  F_fun(state,  dt, out_8208644078465651097);
}
void car_h_25(double *state, double *unused, double *out_4728890029886282683) {
  h_25(state, unused, out_4728890029886282683);
}
void car_H_25(double *state, double *unused, double *out_3148005063319919381) {
  H_25(state, unused, out_3148005063319919381);
}
void car_h_24(double *state, double *unused, double *out_6821372846868724365) {
  h_24(state, unused, out_6821372846868724365);
}
void car_H_24(double *state, double *unused, double *out_3993001672103614780) {
  H_24(state, unused, out_3993001672103614780);
}
void car_h_30(double *state, double *unused, double *out_2090612130942494452) {
  h_30(state, unused, out_2090612130942494452);
}
void car_H_30(double *state, double *unused, double *out_8382048668898015480) {
  H_30(state, unused, out_8382048668898015480);
}
void car_h_26(double *state, double *unused, double *out_233144791967815326) {
  h_26(state, unused, out_233144791967815326);
}
void car_H_26(double *state, double *unused, double *out_6452531033080719982) {
  H_26(state, unused, out_6452531033080719982);
}
void car_h_27(double *state, double *unused, double *out_7850890224567832796) {
  h_27(state, unused, out_7850890224567832796);
}
void car_H_27(double *state, double *unused, double *out_7889932093011111225) {
  H_27(state, unused, out_7889932093011111225);
}
void car_h_29(double *state, double *unused, double *out_1587810802521718556) {
  h_29(state, unused, out_1587810802521718556);
}
void car_H_29(double *state, double *unused, double *out_7871817324583623296) {
  H_29(state, unused, out_7871817324583623296);
}
void car_h_28(double *state, double *unused, double *out_7991448528333021169) {
  h_28(state, unused, out_7991448528333021169);
}
void car_H_28(double *state, double *unused, double *out_1094170349072029618) {
  H_28(state, unused, out_1094170349072029618);
}
void car_h_31(double *state, double *unused, double *out_1732497105052371523) {
  h_31(state, unused, out_1732497105052371523);
}
void car_H_31(double *state, double *unused, double *out_3178651025196879809) {
  H_31(state, unused, out_3178651025196879809);
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
