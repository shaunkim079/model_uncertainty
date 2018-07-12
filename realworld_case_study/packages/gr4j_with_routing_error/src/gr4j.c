// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>


#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b)?(a):(b))
#define max(a, b) ((a < b)?(b):(a))
#endif


void sma_gr4j_sk(double *P, double *E, int *n,
	 double *x1, double *param2, double *initial_state, double *state_error,
         double *U, double *S, double *ET, int *success)
{
    int t;
    double Pn, En, Ps, St_x1, perc;
    double S_prev = *initial_state;
    for (t = 0; t < *n; t++) {
		if (t > 0) {
			S_prev += state_error[t-1];
		}
        Pn = max(P[t] - E[t], 0);
        En = max(E[t] - P[t], 0);
        St_x1 = S_prev / (*x1);
        // production
        Ps = 0;
        ET[t] = 0;
        if (Pn > 0) {
            // part of Pn fills the production store
            Ps = ( *x1 * (1.0 - pow(St_x1, *param2)) * tanh(Pn / (*x1)) /
                   (1.0 + St_x1 * tanh(Pn / (*x1))) );
        }
        if (En > 0) {
            // actual evaporation
            ET[t] = ( S_prev * (2.0 - St_x1) * tanh(En / (*x1)) /
                      (1.0 + (1.0 - St_x1) * tanh(En / (*x1))) );
        }
        S[t] = S_prev - ET[t] + Ps;
        // percolation leakage
        perc = S[t] * ( 1.0 - pow(1.0 + pow((4.0/9.0) * S[t] / (*x1), 4.0), -0.25) );
        S[t] = S[t] - perc;
        U[t] = perc + (Pn - Ps);
        S_prev = S[t];
		
		if (S_prev<0 || U[t]<0 || ET[t]<0 || Ps<0 || perc<0) {
			*success = 0;
			break;
		}
    }
}

void routing_gr4j_sk(double *Q9, double *Q1, int *n,
             double *x2, double *x3, double *R_0,
             double *Qr, double *Qd, double *R, double *state_error_R,
			 int *success)
{
    int t;
    double Rt_x3, F;
    double R_prev = (*R_0) * (*x3);
    for (t = 0; t < *n; t++) {
		if (t > 0) {
			R_prev += state_error_R[t-1];
		}
        Rt_x3 = R_prev / (*x3);
        // groundwater exchange term
        F = (*x2) * pow(Rt_x3, 7.0/2.0);
        // reservoir level
        R[t] = max(0, R_prev + Q9[t] + F);
        // outflow of reservoir
        Qr[t] = R[t] * ( 1.0 - pow(1.0 + pow(R[t] / (*x3), 4.0), -0.25) );
        R[t] = R[t] - Qr[t];
        // other store
        Qd[t] = max(0, Q1[t] + F);
        R_prev = R[t];
		
		if (R_prev<0 || Rt_x3<0) {
			*success = 0;
			break;
		}
    }
}
