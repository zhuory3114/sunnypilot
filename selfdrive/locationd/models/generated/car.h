#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_6387681849755696540);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1447910342006919488);
void car_H_mod_fun(double *state, double *out_3269785016151077926);
void car_f_fun(double *state, double dt, double *out_6160736893164847577);
void car_F_fun(double *state, double dt, double *out_8208644078465651097);
void car_h_25(double *state, double *unused, double *out_4728890029886282683);
void car_H_25(double *state, double *unused, double *out_3148005063319919381);
void car_h_24(double *state, double *unused, double *out_6821372846868724365);
void car_H_24(double *state, double *unused, double *out_3993001672103614780);
void car_h_30(double *state, double *unused, double *out_2090612130942494452);
void car_H_30(double *state, double *unused, double *out_8382048668898015480);
void car_h_26(double *state, double *unused, double *out_233144791967815326);
void car_H_26(double *state, double *unused, double *out_6452531033080719982);
void car_h_27(double *state, double *unused, double *out_7850890224567832796);
void car_H_27(double *state, double *unused, double *out_7889932093011111225);
void car_h_29(double *state, double *unused, double *out_1587810802521718556);
void car_H_29(double *state, double *unused, double *out_7871817324583623296);
void car_h_28(double *state, double *unused, double *out_7991448528333021169);
void car_H_28(double *state, double *unused, double *out_1094170349072029618);
void car_h_31(double *state, double *unused, double *out_1732497105052371523);
void car_H_31(double *state, double *unused, double *out_3178651025196879809);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}