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

/*  Migrate data into migdata */
void mig_kirch_loop (double * data, double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, bool nearfield, double * dist2, double * rs){
    int i;
    int j;
    int k;
    int l;
    double min;
    double m;
    double costheta;
    double integral;
    int Didx;
    int lmin_prev;
    int rmin_prev;
    setbuf(stdout, NULL);
    printf("Iter: ");

    
    for(j=0;j<tnum;j++){
        printf("%d,", j);
        for(i=0;i<snum;i++){
            integral = 0.0;
            rmin_prev = i;
            lmin_prev = i;
            /* Do left and right separately so we know when to break */
            for(k=j + 1;k<tnum;k++){
                rs[k] = sqrt(pow(dist[k] - dist[j], 2.) + zs2[j]);
                if(2. * rs[k] / vel > max_travel_time){
                     break;
                }
                costheta = zs[j] / rs[k];
                min = 1.0e6;
                /* this should not be smaller than the previous rmin? */
                for(l=rmin_prev;l<snum;l++){
                    m = fabs(tt_sec[l] - 2. * rs[l] / vel);
                    if(m < min){
                        min = m;
                        Didx = l;
                    }
                }
                rmin_prev = l;
                integral += gradD[k * snum + Didx] * costheta / vel;
            }

            for(k=j;k>=0;k--){
                rs[k] = sqrt(pow(dist[k] - dist[j], 2.) + zs2[j]);
                if(2. * rs[k] / vel > max_travel_time){
                     break;
                }
                costheta = zs[j] / rs[k];
                min = 1.0e6;
                /* this should not be smaller than the previous lmin? */
                for(l=lmin_prev;l<snum;l++){
                    m = fabs(tt_sec[l] - 2. * rs[l] / vel);
                    if(m < min){
                        min = m;
                        Didx = l;
                    }
                }
                lmin_prev = l;
                integral += gradD[k * snum + Didx] * costheta / vel;
            }
            migdata[j * snum + i] = integral;
        }
    }
    printf("\n");
}
