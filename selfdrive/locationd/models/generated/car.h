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
void car_err_fun(double *nom_x, double *delta_x, double *out_8677092309422431987);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3536033552698117999);
void car_H_mod_fun(double *state, double *out_8253580055064896978);
void car_f_fun(double *state, double dt, double *out_5986056227021830781);
void car_F_fun(double *state, double dt, double *out_8361123820946526770);
void car_h_25(double *state, double *unused, double *out_7030946204566009936);
void car_H_25(double *state, double *unused, double *out_3735175456469797709);
void car_h_24(double *state, double *unused, double *out_1397172579223915339);
void car_H_24(double *state, double *unused, double *out_8900182986216976942);
void car_h_30(double *state, double *unused, double *out_3393496676311987586);
void car_H_30(double *state, double *unused, double *out_7794878275748137152);
void car_h_26(double *state, double *unused, double *out_5588472495694853179);
void car_H_26(double *state, double *unused, double *out_7039701426230598310);
void car_h_27(double *state, double *unused, double *out_2029272936905164291);
void car_H_27(double *state, double *unused, double *out_8477102486160989553);
void car_h_29(double *state, double *unused, double *out_1608176591348858059);
void car_H_29(double *state, double *unused, double *out_7284646931433744968);
void car_h_28(double *state, double *unused, double *out_5669780174728713609);
void car_H_28(double *state, double *unused, double *out_6079698125206276074);
void car_h_31(double *state, double *unused, double *out_7554114392715040028);
void car_H_31(double *state, double *unused, double *out_3765821418346758137);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}