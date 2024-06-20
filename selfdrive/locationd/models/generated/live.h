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
void live_H(double *in_vec, double *out_1247131383871684740);
void live_err_fun(double *nom_x, double *delta_x, double *out_4876428447537604078);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6213060246611821173);
void live_H_mod_fun(double *state, double *out_8355549558907757707);
void live_f_fun(double *state, double dt, double *out_4469768777619126135);
void live_F_fun(double *state, double dt, double *out_3099738278590939435);
void live_h_4(double *state, double *unused, double *out_4018366025072394088);
void live_H_4(double *state, double *unused, double *out_6578178129224055764);
void live_h_9(double *state, double *unused, double *out_7591538521268605337);
void live_H_9(double *state, double *unused, double *out_6336988482594465119);
void live_h_10(double *state, double *unused, double *out_6753538431309370865);
void live_H_10(double *state, double *unused, double *out_433116642984162329);
void live_h_12(double *state, double *unused, double *out_439037189819194306);
void live_H_12(double *state, double *unused, double *out_1558721721192093969);
void live_h_35(double *state, double *unused, double *out_6513537197017262713);
void live_H_35(double *state, double *unused, double *out_3211516071851448388);
void live_h_32(double *state, double *unused, double *out_5491330579023163217);
void live_H_32(double *state, double *unused, double *out_6097997494699902884);
void live_h_13(double *state, double *unused, double *out_3444580219625736867);
void live_H_13(double *state, double *unused, double *out_1852456986602090872);
void live_h_14(double *state, double *unused, double *out_7591538521268605337);
void live_H_14(double *state, double *unused, double *out_6336988482594465119);
void live_h_33(double *state, double *unused, double *out_1713139368494136663);
void live_H_33(double *state, double *unused, double *out_60959067212590784);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}