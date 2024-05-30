#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_8848504734180065506);
void live_err_fun(double *nom_x, double *delta_x, double *out_1831232880786325784);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1821126222288460850);
void live_H_mod_fun(double *state, double *out_3908414896028627082);
void live_f_fun(double *state, double dt, double *out_7892289561224452355);
void live_F_fun(double *state, double dt, double *out_8515547901378915655);
void live_h_4(double *state, double *unused, double *out_8679006882303680840);
void live_H_4(double *state, double *unused, double *out_5990756093177091639);
void live_h_9(double *state, double *unused, double *out_2409868858877289513);
void live_H_9(double *state, double *unused, double *out_1296462842087355831);
void live_h_10(double *state, double *unused, double *out_120384295002916867);
void live_H_10(double *state, double *unused, double *out_3553941899114977873);
void live_h_12(double *state, double *unused, double *out_7598365100289183419);
void live_H_12(double *state, double *unused, double *out_1676372220505358853);
void live_h_35(double *state, double *unused, double *out_2019650774056222022);
void live_H_35(double *state, double *unused, double *out_4421935252830372562);
void live_h_32(double *state, double *unused, double *out_5076051470792844394);
void live_H_32(double *state, double *unused, double *out_2492929250863743794);
void live_h_13(double *state, double *unused, double *out_3672272796614809703);
void live_H_13(double *state, double *unused, double *out_8928872552718057093);
void live_h_14(double *state, double *unused, double *out_2409868858877289513);
void live_H_14(double *state, double *unused, double *out_1296462842087355831);
void live_h_33(double *state, double *unused, double *out_9121223161723803206);
void live_H_33(double *state, double *unused, double *out_7572492257469230166);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}