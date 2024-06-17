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
void live_H(double *in_vec, double *out_775903400974484091);
void live_err_fun(double *nom_x, double *delta_x, double *out_4354435775177582705);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_7121017768917172192);
void live_H_mod_fun(double *state, double *out_8943981757272317000);
void live_f_fun(double *state, double dt, double *out_2068157913889425259);
void live_F_fun(double *state, double dt, double *out_7057504128973734560);
void live_h_4(double *state, double *unused, double *out_1466244038284864942);
void live_H_4(double *state, double *unused, double *out_4719075018422975890);
void live_h_9(double *state, double *unused, double *out_4645040410163205132);
void live_H_9(double *state, double *unused, double *out_4960264665052566535);
void live_h_10(double *state, double *unused, double *out_4702588874207177269);
void live_H_10(double *state, double *unused, double *out_5625499360276871869);
void live_h_12(double *state, double *unused, double *out_3263115080229759064);
void live_H_12(double *state, double *unused, double *out_8708212647254613931);
void live_h_35(double *state, double *unused, double *out_2379059252865934280);
void live_H_35(double *state, double *unused, double *out_8085737075795583266);
void live_h_32(double *state, double *unused, double *out_2251998733377667084);
void live_H_32(double *state, double *unused, double *out_4559702625706355792);
void live_h_13(double *state, double *unused, double *out_5340354336633023650);
void live_H_13(double *state, double *unused, double *out_7408247029735667991);
void live_h_14(double *state, double *unused, double *out_4645040410163205132);
void live_H_14(double *state, double *unused, double *out_4960264665052566535);
void live_h_33(double *state, double *unused, double *out_4121750344094741723);
void live_H_33(double *state, double *unused, double *out_7210449993275110746);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}