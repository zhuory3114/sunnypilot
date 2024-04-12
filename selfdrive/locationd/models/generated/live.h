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
void live_H(double *in_vec, double *out_7354899081148072319);
void live_err_fun(double *nom_x, double *delta_x, double *out_926748312842793590);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8553861514595348833);
void live_H_mod_fun(double *state, double *out_1090406756321095952);
void live_f_fun(double *state, double dt, double *out_1607575285317213045);
void live_F_fun(double *state, double dt, double *out_3881758448010921799);
void live_h_4(double *state, double *unused, double *out_5351864317507456224);
void live_H_4(double *state, double *unused, double *out_9128846745444311004);
void live_h_9(double *state, double *unused, double *out_2329374119852845347);
void live_H_9(double *state, double *unused, double *out_2030678393000793142);
void live_h_10(double *state, double *unused, double *out_2828269157780052031);
void live_H_10(double *state, double *unused, double *out_1320877861754912006);
void live_h_12(double *state, double *unused, double *out_6763182500966395644);
void live_H_12(double *state, double *unused, double *out_1650769014582790120);
void live_h_35(double *state, double *unused, double *out_7964897784483847122);
void live_H_35(double *state, double *unused, double *out_5951235270892633236);
void live_h_32(double *state, double *unused, double *out_415144153704329373);
void live_H_32(double *state, double *unused, double *out_4642638357729382318);
void live_h_13(double *state, double *unused, double *out_8033899621259556965);
void live_H_13(double *state, double *unused, double *out_1251369415488548331);
void live_h_14(double *state, double *unused, double *out_2329374119852845347);
void live_H_14(double *state, double *unused, double *out_2030678393000793142);
void live_h_33(double *state, double *unused, double *out_7785328566767514480);
void live_H_33(double *state, double *unused, double *out_2800678266253775632);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}