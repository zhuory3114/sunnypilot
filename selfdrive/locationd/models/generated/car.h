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
void car_err_fun(double *nom_x, double *delta_x, double *out_3930116420352682607);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8332055307244711483);
void car_H_mod_fun(double *state, double *out_1183099394398996416);
void car_f_fun(double *state, double dt, double *out_3974371446758158566);
void car_F_fun(double *state, double dt, double *out_1509814126693874228);
void car_h_25(double *state, double *unused, double *out_479281131967255891);
void car_H_25(double *state, double *unused, double *out_2945433631051225152);
void car_h_24(double *state, double *unused, double *out_2798532925715571291);
void car_H_24(double *state, double *unused, double *out_772784032045725586);
void car_h_30(double *state, double *unused, double *out_2517111792866655269);
void car_H_30(double *state, double *unused, double *out_5463766589558473779);
void car_h_26(double *state, double *unused, double *out_6194954447731099144);
void car_H_26(double *state, double *unused, double *out_796069687822831072);
void car_h_27(double *state, double *unused, double *out_1422790984868918326);
void car_H_27(double *state, double *unused, double *out_3757026010876807957);
void car_h_29(double *state, double *unused, double *out_8159841663947612104);
void car_H_29(double *state, double *unused, double *out_1072031354761990862);
void car_h_28(double *state, double *unused, double *out_5704354826112578834);
void car_H_28(double *state, double *unused, double *out_891598916803335389);
void car_h_31(double *state, double *unused, double *out_8360553782721305076);
void car_H_31(double *state, double *unused, double *out_1422277790056182548);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}