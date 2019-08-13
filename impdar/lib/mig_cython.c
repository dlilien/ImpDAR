/*
 * mig_cython.c
 * Copyright (C) 2019 dlilien <dlilien@berens>
 *
 * Distributed under terms of the MIT license.
 */

#include "mig_cython.h"
#include <stdbool.h>


#include <math.h>
#include <stdio.h>

/*  Compute the cosine of each element in in_array, storing the result in
 *  out_array. */
void mig_cython (double * data, double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, bool nearfield){
    int i;
    int j;
    int k;
    int l;
    double min;
    double m;
    double dist2[tnum];
    double rs[snum];
    double costheta;
    double integral;
    int Didx;
    double gradDhyp;
    for(j=0;j<tnum;j++){
        for(i=0;i<tnum;i++){
            dist2[i] = pow(dist[i] - dist[j], 2.);
        }
        for(i=0;i<snum;i++){
            for(k=0;k<snum;k++){
                rs[k] = sqrt(dist2[k] + zs2[j]);
            }
            printf("rs is calculated at iter: %d\n", i);
            integral = 0.0;
            for(k=0;k<tnum;k++){
                printf("Inside iter: %d\n", k);
                costheta = zs[j] / rs[k];
                min = 1.0e6;
                for(l=0;l<snum;l++){
                    printf("Inside iter: %d\n", l);
                    m = tt_sec[l] - 2. * rs[l] / vel;
                    printf("%f %f %d\n", m, min, Didx);
                    if(m < min){
                        min = m;
                        Didx = l;
                    }
                }
                gradDhyp = gradD[Didx * snum + k];
                if(2. * rs[k] / vel > max_travel_time){
                    gradDhyp=0.;
                }
                integral += gradDhyp * costheta / vel;
                printf("Finished inside loop\n");
            }
            migdata[j * snum + i] = integral;
        }
    }
}
