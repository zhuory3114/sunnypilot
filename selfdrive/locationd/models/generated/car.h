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
void car_err_fun(double *nom_x, double *delta_x, double *out_2907344304532562836);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8954218720378571779);
void car_H_mod_fun(double *state, double *out_6798369121237450864);
void car_f_fun(double *state, double dt, double *out_6734862407766948837);
void car_F_fun(double *state, double dt, double *out_7156118317797004448);
void car_h_25(double *state, double *unused, double *out_5350053244966229021);
void car_H_25(double *state, double *unused, double *out_2711664023153060623);
void car_h_24(double *state, double *unused, double *out_5232067606711025854);
void car_H_24(double *state, double *unused, double *out_539014424147561057);
void car_h_30(double *state, double *unused, double *out_5636132064494434270);
void car_H_30(double *state, double *unused, double *out_5229996981660309250);
void car_h_26(double *state, double *unused, double *out_6554787287256070092);
void car_H_26(double *state, double *unused, double *out_1029839295720995601);
void car_h_27(double *state, double *unused, double *out_3670490167474247421);
void car_H_27(double *state, double *unused, double *out_3055233669859884339);
void car_h_29(double *state, double *unused, double *out_3513303499123118276);
void car_H_29(double *state, double *unused, double *out_5740228325974701434);
void car_h_28(double *state, double *unused, double *out_2728551349939601725);
void car_H_28(double *state, double *unused, double *out_657829308905170860);
void car_h_31(double *state, double *unused, double *out_8987502773220251371);
void car_H_31(double *state, double *unused, double *out_1656047397954347077);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}