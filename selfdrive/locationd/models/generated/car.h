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
void car_err_fun(double *nom_x, double *delta_x, double *out_6568424505483499173);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4235031182312661852);
void car_H_mod_fun(double *state, double *out_1527419437699775404);
void car_f_fun(double *state, double dt, double *out_2784621583213052968);
void car_F_fun(double *state, double dt, double *out_2545055737197701349);
void car_h_25(double *state, double *unused, double *out_8039810018730752828);
void car_H_25(double *state, double *unused, double *out_1131102053972815576);
void car_h_24(double *state, double *unused, double *out_7458370362274233029);
void car_H_24(double *state, double *unused, double *out_5877611380939324610);
void car_h_30(double *state, double *unused, double *out_1599795096106713024);
void car_H_30(double *state, double *unused, double *out_5785588287518801179);
void car_h_26(double *state, double *unused, double *out_8441292157958611939);
void car_H_26(double *state, double *unused, double *out_4872605372846871800);
void car_h_27(double *state, double *unused, double *out_3077430021049459596);
void car_H_27(double *state, double *unused, double *out_3610824975718376268);
void car_h_29(double *state, double *unused, double *out_3885288935234119431);
void car_H_29(double *state, double *unused, double *out_6295819631833193363);
void car_h_28(double *state, double *unused, double *out_473024700240761945);
void car_H_28(double *state, double *unused, double *out_1213420614763662789);
void car_h_31(double *state, double *unused, double *out_2631300804724111224);
void car_H_31(double *state, double *unused, double *out_1100456092095855148);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}