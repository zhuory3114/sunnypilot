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
void err_fun(double *nom_x, double *delta_x, double *out_6568424505483499173) {
   out_6568424505483499173[0] = delta_x[0] + nom_x[0];
   out_6568424505483499173[1] = delta_x[1] + nom_x[1];
   out_6568424505483499173[2] = delta_x[2] + nom_x[2];
   out_6568424505483499173[3] = delta_x[3] + nom_x[3];
   out_6568424505483499173[4] = delta_x[4] + nom_x[4];
   out_6568424505483499173[5] = delta_x[5] + nom_x[5];
   out_6568424505483499173[6] = delta_x[6] + nom_x[6];
   out_6568424505483499173[7] = delta_x[7] + nom_x[7];
   out_6568424505483499173[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4235031182312661852) {
   out_4235031182312661852[0] = -nom_x[0] + true_x[0];
   out_4235031182312661852[1] = -nom_x[1] + true_x[1];
   out_4235031182312661852[2] = -nom_x[2] + true_x[2];
   out_4235031182312661852[3] = -nom_x[3] + true_x[3];
   out_4235031182312661852[4] = -nom_x[4] + true_x[4];
   out_4235031182312661852[5] = -nom_x[5] + true_x[5];
   out_4235031182312661852[6] = -nom_x[6] + true_x[6];
   out_4235031182312661852[7] = -nom_x[7] + true_x[7];
   out_4235031182312661852[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1527419437699775404) {
   out_1527419437699775404[0] = 1.0;
   out_1527419437699775404[1] = 0;
   out_1527419437699775404[2] = 0;
   out_1527419437699775404[3] = 0;
   out_1527419437699775404[4] = 0;
   out_1527419437699775404[5] = 0;
   out_1527419437699775404[6] = 0;
   out_1527419437699775404[7] = 0;
   out_1527419437699775404[8] = 0;
   out_1527419437699775404[9] = 0;
   out_1527419437699775404[10] = 1.0;
   out_1527419437699775404[11] = 0;
   out_1527419437699775404[12] = 0;
   out_1527419437699775404[13] = 0;
   out_1527419437699775404[14] = 0;
   out_1527419437699775404[15] = 0;
   out_1527419437699775404[16] = 0;
   out_1527419437699775404[17] = 0;
   out_1527419437699775404[18] = 0;
   out_1527419437699775404[19] = 0;
   out_1527419437699775404[20] = 1.0;
   out_1527419437699775404[21] = 0;
   out_1527419437699775404[22] = 0;
   out_1527419437699775404[23] = 0;
   out_1527419437699775404[24] = 0;
   out_1527419437699775404[25] = 0;
   out_1527419437699775404[26] = 0;
   out_1527419437699775404[27] = 0;
   out_1527419437699775404[28] = 0;
   out_1527419437699775404[29] = 0;
   out_1527419437699775404[30] = 1.0;
   out_1527419437699775404[31] = 0;
   out_1527419437699775404[32] = 0;
   out_1527419437699775404[33] = 0;
   out_1527419437699775404[34] = 0;
   out_1527419437699775404[35] = 0;
   out_1527419437699775404[36] = 0;
   out_1527419437699775404[37] = 0;
   out_1527419437699775404[38] = 0;
   out_1527419437699775404[39] = 0;
   out_1527419437699775404[40] = 1.0;
   out_1527419437699775404[41] = 0;
   out_1527419437699775404[42] = 0;
   out_1527419437699775404[43] = 0;
   out_1527419437699775404[44] = 0;
   out_1527419437699775404[45] = 0;
   out_1527419437699775404[46] = 0;
   out_1527419437699775404[47] = 0;
   out_1527419437699775404[48] = 0;
   out_1527419437699775404[49] = 0;
   out_1527419437699775404[50] = 1.0;
   out_1527419437699775404[51] = 0;
   out_1527419437699775404[52] = 0;
   out_1527419437699775404[53] = 0;
   out_1527419437699775404[54] = 0;
   out_1527419437699775404[55] = 0;
   out_1527419437699775404[56] = 0;
   out_1527419437699775404[57] = 0;
   out_1527419437699775404[58] = 0;
   out_1527419437699775404[59] = 0;
   out_1527419437699775404[60] = 1.0;
   out_1527419437699775404[61] = 0;
   out_1527419437699775404[62] = 0;
   out_1527419437699775404[63] = 0;
   out_1527419437699775404[64] = 0;
   out_1527419437699775404[65] = 0;
   out_1527419437699775404[66] = 0;
   out_1527419437699775404[67] = 0;
   out_1527419437699775404[68] = 0;
   out_1527419437699775404[69] = 0;
   out_1527419437699775404[70] = 1.0;
   out_1527419437699775404[71] = 0;
   out_1527419437699775404[72] = 0;
   out_1527419437699775404[73] = 0;
   out_1527419437699775404[74] = 0;
   out_1527419437699775404[75] = 0;
   out_1527419437699775404[76] = 0;
   out_1527419437699775404[77] = 0;
   out_1527419437699775404[78] = 0;
   out_1527419437699775404[79] = 0;
   out_1527419437699775404[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2784621583213052968) {
   out_2784621583213052968[0] = state[0];
   out_2784621583213052968[1] = state[1];
   out_2784621583213052968[2] = state[2];
   out_2784621583213052968[3] = state[3];
   out_2784621583213052968[4] = state[4];
   out_2784621583213052968[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2784621583213052968[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2784621583213052968[7] = state[7];
   out_2784621583213052968[8] = state[8];
}
void F_fun(double *state, double dt, double *out_2545055737197701349) {
   out_2545055737197701349[0] = 1;
   out_2545055737197701349[1] = 0;
   out_2545055737197701349[2] = 0;
   out_2545055737197701349[3] = 0;
   out_2545055737197701349[4] = 0;
   out_2545055737197701349[5] = 0;
   out_2545055737197701349[6] = 0;
   out_2545055737197701349[7] = 0;
   out_2545055737197701349[8] = 0;
   out_2545055737197701349[9] = 0;
   out_2545055737197701349[10] = 1;
   out_2545055737197701349[11] = 0;
   out_2545055737197701349[12] = 0;
   out_2545055737197701349[13] = 0;
   out_2545055737197701349[14] = 0;
   out_2545055737197701349[15] = 0;
   out_2545055737197701349[16] = 0;
   out_2545055737197701349[17] = 0;
   out_2545055737197701349[18] = 0;
   out_2545055737197701349[19] = 0;
   out_2545055737197701349[20] = 1;
   out_2545055737197701349[21] = 0;
   out_2545055737197701349[22] = 0;
   out_2545055737197701349[23] = 0;
   out_2545055737197701349[24] = 0;
   out_2545055737197701349[25] = 0;
   out_2545055737197701349[26] = 0;
   out_2545055737197701349[27] = 0;
   out_2545055737197701349[28] = 0;
   out_2545055737197701349[29] = 0;
   out_2545055737197701349[30] = 1;
   out_2545055737197701349[31] = 0;
   out_2545055737197701349[32] = 0;
   out_2545055737197701349[33] = 0;
   out_2545055737197701349[34] = 0;
   out_2545055737197701349[35] = 0;
   out_2545055737197701349[36] = 0;
   out_2545055737197701349[37] = 0;
   out_2545055737197701349[38] = 0;
   out_2545055737197701349[39] = 0;
   out_2545055737197701349[40] = 1;
   out_2545055737197701349[41] = 0;
   out_2545055737197701349[42] = 0;
   out_2545055737197701349[43] = 0;
   out_2545055737197701349[44] = 0;
   out_2545055737197701349[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2545055737197701349[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2545055737197701349[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2545055737197701349[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2545055737197701349[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2545055737197701349[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2545055737197701349[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2545055737197701349[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2545055737197701349[53] = -9.8000000000000007*dt;
   out_2545055737197701349[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2545055737197701349[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2545055737197701349[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2545055737197701349[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2545055737197701349[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2545055737197701349[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2545055737197701349[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2545055737197701349[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2545055737197701349[62] = 0;
   out_2545055737197701349[63] = 0;
   out_2545055737197701349[64] = 0;
   out_2545055737197701349[65] = 0;
   out_2545055737197701349[66] = 0;
   out_2545055737197701349[67] = 0;
   out_2545055737197701349[68] = 0;
   out_2545055737197701349[69] = 0;
   out_2545055737197701349[70] = 1;
   out_2545055737197701349[71] = 0;
   out_2545055737197701349[72] = 0;
   out_2545055737197701349[73] = 0;
   out_2545055737197701349[74] = 0;
   out_2545055737197701349[75] = 0;
   out_2545055737197701349[76] = 0;
   out_2545055737197701349[77] = 0;
   out_2545055737197701349[78] = 0;
   out_2545055737197701349[79] = 0;
   out_2545055737197701349[80] = 1;
}
void h_25(double *state, double *unused, double *out_8039810018730752828) {
   out_8039810018730752828[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1131102053972815576) {
   out_1131102053972815576[0] = 0;
   out_1131102053972815576[1] = 0;
   out_1131102053972815576[2] = 0;
   out_1131102053972815576[3] = 0;
   out_1131102053972815576[4] = 0;
   out_1131102053972815576[5] = 0;
   out_1131102053972815576[6] = 1;
   out_1131102053972815576[7] = 0;
   out_1131102053972815576[8] = 0;
}
void h_24(double *state, double *unused, double *out_7458370362274233029) {
   out_7458370362274233029[0] = state[4];
   out_7458370362274233029[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5877611380939324610) {
   out_5877611380939324610[0] = 0;
   out_5877611380939324610[1] = 0;
   out_5877611380939324610[2] = 0;
   out_5877611380939324610[3] = 0;
   out_5877611380939324610[4] = 1;
   out_5877611380939324610[5] = 0;
   out_5877611380939324610[6] = 0;
   out_5877611380939324610[7] = 0;
   out_5877611380939324610[8] = 0;
   out_5877611380939324610[9] = 0;
   out_5877611380939324610[10] = 0;
   out_5877611380939324610[11] = 0;
   out_5877611380939324610[12] = 0;
   out_5877611380939324610[13] = 0;
   out_5877611380939324610[14] = 1;
   out_5877611380939324610[15] = 0;
   out_5877611380939324610[16] = 0;
   out_5877611380939324610[17] = 0;
}
void h_30(double *state, double *unused, double *out_1599795096106713024) {
   out_1599795096106713024[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5785588287518801179) {
   out_5785588287518801179[0] = 0;
   out_5785588287518801179[1] = 0;
   out_5785588287518801179[2] = 0;
   out_5785588287518801179[3] = 0;
   out_5785588287518801179[4] = 1;
   out_5785588287518801179[5] = 0;
   out_5785588287518801179[6] = 0;
   out_5785588287518801179[7] = 0;
   out_5785588287518801179[8] = 0;
}
void h_26(double *state, double *unused, double *out_8441292157958611939) {
   out_8441292157958611939[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4872605372846871800) {
   out_4872605372846871800[0] = 0;
   out_4872605372846871800[1] = 0;
   out_4872605372846871800[2] = 0;
   out_4872605372846871800[3] = 0;
   out_4872605372846871800[4] = 0;
   out_4872605372846871800[5] = 0;
   out_4872605372846871800[6] = 0;
   out_4872605372846871800[7] = 1;
   out_4872605372846871800[8] = 0;
}
void h_27(double *state, double *unused, double *out_3077430021049459596) {
   out_3077430021049459596[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3610824975718376268) {
   out_3610824975718376268[0] = 0;
   out_3610824975718376268[1] = 0;
   out_3610824975718376268[2] = 0;
   out_3610824975718376268[3] = 1;
   out_3610824975718376268[4] = 0;
   out_3610824975718376268[5] = 0;
   out_3610824975718376268[6] = 0;
   out_3610824975718376268[7] = 0;
   out_3610824975718376268[8] = 0;
}
void h_29(double *state, double *unused, double *out_3885288935234119431) {
   out_3885288935234119431[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6295819631833193363) {
   out_6295819631833193363[0] = 0;
   out_6295819631833193363[1] = 1;
   out_6295819631833193363[2] = 0;
   out_6295819631833193363[3] = 0;
   out_6295819631833193363[4] = 0;
   out_6295819631833193363[5] = 0;
   out_6295819631833193363[6] = 0;
   out_6295819631833193363[7] = 0;
   out_6295819631833193363[8] = 0;
}
void h_28(double *state, double *unused, double *out_473024700240761945) {
   out_473024700240761945[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1213420614763662789) {
   out_1213420614763662789[0] = 1;
   out_1213420614763662789[1] = 0;
   out_1213420614763662789[2] = 0;
   out_1213420614763662789[3] = 0;
   out_1213420614763662789[4] = 0;
   out_1213420614763662789[5] = 0;
   out_1213420614763662789[6] = 0;
   out_1213420614763662789[7] = 0;
   out_1213420614763662789[8] = 0;
}
void h_31(double *state, double *unused, double *out_2631300804724111224) {
   out_2631300804724111224[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1100456092095855148) {
   out_1100456092095855148[0] = 0;
   out_1100456092095855148[1] = 0;
   out_1100456092095855148[2] = 0;
   out_1100456092095855148[3] = 0;
   out_1100456092095855148[4] = 0;
   out_1100456092095855148[5] = 0;
   out_1100456092095855148[6] = 0;
   out_1100456092095855148[7] = 0;
   out_1100456092095855148[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_6568424505483499173) {
  err_fun(nom_x, delta_x, out_6568424505483499173);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4235031182312661852) {
  inv_err_fun(nom_x, true_x, out_4235031182312661852);
}
void car_H_mod_fun(double *state, double *out_1527419437699775404) {
  H_mod_fun(state, out_1527419437699775404);
}
void car_f_fun(double *state, double dt, double *out_2784621583213052968) {
  f_fun(state,  dt, out_2784621583213052968);
}
void car_F_fun(double *state, double dt, double *out_2545055737197701349) {
  F_fun(state,  dt, out_2545055737197701349);
}
void car_h_25(double *state, double *unused, double *out_8039810018730752828) {
  h_25(state, unused, out_8039810018730752828);
}
void car_H_25(double *state, double *unused, double *out_1131102053972815576) {
  H_25(state, unused, out_1131102053972815576);
}
void car_h_24(double *state, double *unused, double *out_7458370362274233029) {
  h_24(state, unused, out_7458370362274233029);
}
void car_H_24(double *state, double *unused, double *out_5877611380939324610) {
  H_24(state, unused, out_5877611380939324610);
}
void car_h_30(double *state, double *unused, double *out_1599795096106713024) {
  h_30(state, unused, out_1599795096106713024);
}
void car_H_30(double *state, double *unused, double *out_5785588287518801179) {
  H_30(state, unused, out_5785588287518801179);
}
void car_h_26(double *state, double *unused, double *out_8441292157958611939) {
  h_26(state, unused, out_8441292157958611939);
}
void car_H_26(double *state, double *unused, double *out_4872605372846871800) {
  H_26(state, unused, out_4872605372846871800);
}
void car_h_27(double *state, double *unused, double *out_3077430021049459596) {
  h_27(state, unused, out_3077430021049459596);
}
void car_H_27(double *state, double *unused, double *out_3610824975718376268) {
  H_27(state, unused, out_3610824975718376268);
}
void car_h_29(double *state, double *unused, double *out_3885288935234119431) {
  h_29(state, unused, out_3885288935234119431);
}
void car_H_29(double *state, double *unused, double *out_6295819631833193363) {
  H_29(state, unused, out_6295819631833193363);
}
void car_h_28(double *state, double *unused, double *out_473024700240761945) {
  h_28(state, unused, out_473024700240761945);
}
void car_H_28(double *state, double *unused, double *out_1213420614763662789) {
  H_28(state, unused, out_1213420614763662789);
}
void car_h_31(double *state, double *unused, double *out_2631300804724111224) {
  h_31(state, unused, out_2631300804724111224);
}
void car_H_31(double *state, double *unused, double *out_1100456092095855148) {
  H_31(state, unused, out_1100456092095855148);
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
