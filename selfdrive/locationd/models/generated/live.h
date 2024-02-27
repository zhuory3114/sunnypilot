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
void live_H(double *in_vec, double *out_3793607795253229226);
void live_err_fun(double *nom_x, double *delta_x, double *out_7452277575920999010);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_9084836373951615154);
void live_H_mod_fun(double *state, double *out_3476564020925666246);
void live_f_fun(double *state, double dt, double *out_8632433294484326562);
void live_F_fun(double *state, double dt, double *out_755428289952796010);
void live_h_4(double *state, double *unused, double *out_4135544437234870366);
void live_H_4(double *state, double *unused, double *out_6182807741270379125);
void live_h_9(double *state, double *unused, double *out_40914151002627374);
void live_H_9(double *state, double *unused, double *out_9071669293550458467);
void live_h_10(double *state, double *unused, double *out_6389008795956046356);
void live_H_10(double *state, double *unused, double *out_472665177758661559);
void live_h_12(double *state, double *unused, double *out_7639163368799458896);
void live_H_12(double *state, double *unused, double *out_4596808018756721999);
void live_h_35(double *state, double *unused, double *out_1234611402611119700);
void live_H_35(double *state, double *unused, double *out_8897274275066565115);
void live_h_32(double *state, double *unused, double *out_1720010960892601322);
void live_H_32(double *state, double *unused, double *out_932413177986709303);
void live_h_13(double *state, double *unused, double *out_9098570744165599777);
void live_H_13(double *state, double *unused, double *out_5445158903823703855);
void live_h_14(double *state, double *unused, double *out_40914151002627374);
void live_H_14(double *state, double *unused, double *out_9071669293550458467);
void live_h_33(double *state, double *unused, double *out_8139287745808185532);
void live_H_33(double *state, double *unused, double *out_5746717270427707511);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}