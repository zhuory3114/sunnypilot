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
void live_H(double *in_vec, double *out_1111382426087068543);
void live_err_fun(double *nom_x, double *delta_x, double *out_1764225867510345998);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7174499109594114451);
void live_H_mod_fun(double *state, double *out_224910195578200458);
void live_f_fun(double *state, double dt, double *out_5819282905458513110);
void live_F_fun(double *state, double dt, double *out_3037666247417097889);
void live_h_4(double *state, double *unused, double *out_7244893405169961555);
void live_H_4(double *state, double *unused, double *out_3786420766390017315);
void live_h_9(double *state, double *unused, double *out_4612933811799353023);
void live_H_9(double *state, double *unused, double *out_7373104372055086831);
void live_h_10(double *state, double *unused, double *out_281047318752624896);
void live_H_10(double *state, double *unused, double *out_7713613235007542673);
void live_h_12(double *state, double *unused, double *out_8430576866625388783);
void live_H_12(double *state, double *unused, double *out_2594837610652715681);
void live_h_35(double *state, double *unused, double *out_2164896017631207028);
void live_H_35(double *state, double *unused, double *out_4247631961312070100);
void live_h_32(double *state, double *unused, double *out_3911137531967404423);
void live_H_32(double *state, double *unused, double *out_1553672410895542458);
void live_h_13(double *state, double *unused, double *out_8272517388556683793);
void live_H_13(double *state, double *unused, double *out_9108738720012350504);
void live_h_14(double *state, double *unused, double *out_4612933811799353023);
void live_H_14(double *state, double *unused, double *out_7373104372055086831);
void live_h_33(double *state, double *unused, double *out_5156899043507898016);
void live_H_33(double *state, double *unused, double *out_1097074956673212496);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}