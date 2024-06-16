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
void live_H(double *in_vec, double *out_7093419124625442365);
void live_err_fun(double *nom_x, double *delta_x, double *out_8122962677042080884);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_2879191858436194061);
void live_H_mod_fun(double *state, double *out_6891225504689681335);
void live_f_fun(double *state, double dt, double *out_8358383459853452579);
void live_F_fun(double *state, double dt, double *out_5351070853746756336);
void live_h_4(double *state, double *unused, double *out_1663954682400385106);
void live_H_4(double *state, double *unused, double *out_4643929226346005397);
void live_h_9(double *state, double *unused, double *out_6878549429972817913);
void live_H_9(double *state, double *unused, double *out_6515595912099098749);
void live_h_10(double *state, double *unused, double *out_188068464861658843);
void live_H_10(double *state, double *unused, double *out_4625396864627181921);
void live_h_12(double *state, double *unused, double *out_5577728301545643260);
void live_H_12(double *state, double *unused, double *out_1737329150696727599);
void live_h_35(double *state, double *unused, double *out_352367680913121423);
void live_H_35(double *state, double *unused, double *out_6037795407006570715);
void live_h_32(double *state, double *unused, double *out_6025154312532695622);
void live_H_32(double *state, double *unused, double *out_6350417101894015743);
void live_h_13(double *state, double *unused, double *out_6300134216330282240);
void live_H_13(double *state, double *unused, double *out_735835764300229024);
void live_h_14(double *state, double *unused, double *out_6878549429972817913);
void live_H_14(double *state, double *unused, double *out_6515595912099098749);
void live_h_33(double *state, double *unused, double *out_7772855449669035198);
void live_H_33(double *state, double *unused, double *out_239566496717224414);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}