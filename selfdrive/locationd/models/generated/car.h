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
void car_err_fun(double *nom_x, double *delta_x, double *out_6221613173377875958);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1873114568696790253);
void car_H_mod_fun(double *state, double *out_9120391434276283067);
void car_f_fun(double *state, double dt, double *out_8908813278594628189);
void car_F_fun(double *state, double dt, double *out_285691373668569915);
void car_h_25(double *state, double *unused, double *out_4174992519817443698);
void car_H_25(double *state, double *unused, double *out_2624001866442726016);
void car_h_24(double *state, double *unused, double *out_9139150395481525125);
void car_H_24(double *state, double *unused, double *out_4801216290049875989);
void car_h_30(double *state, double *unused, double *out_3588272737003228743);
void car_H_30(double *state, double *unused, double *out_5142334824949974643);
void car_h_26(double *state, double *unused, double *out_551486685558896058);
void car_H_26(double *state, double *unused, double *out_1117501452431330208);
void car_h_27(double *state, double *unused, double *out_8169232118158913528);
void car_H_27(double *state, double *unused, double *out_2967571513149549732);
void car_h_29(double *state, double *unused, double *out_1413400530657616902);
void car_H_29(double *state, double *unused, double *out_5652566169264366827);
void car_h_28(double *state, double *unused, double *out_7240725904106305403);
void car_H_28(double *state, double *unused, double *out_570167152194836253);
void car_h_31(double *state, double *unused, double *out_5775090872628367926);
void car_H_31(double *state, double *unused, double *out_1743709554664681684);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}