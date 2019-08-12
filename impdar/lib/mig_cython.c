/*
 * mig_cython.c
 * Copyright (C) 2019 dlilien <dlilien@berens>
 *
 * Distributed under terms of the MIT license.
 */

#include "mig_cython.h"
#include <stdbool.h>


#include <math.h>

/*  Compute the cosine of each element in in_array, storing the result in
 *  out_array. */
void mig_cython (double * data, double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, bool nearfield){
    int i;
    int j;
    j = 0;
    for(i=0;i<snum;i++){
        migdata[i] = cos(data[i]);
    }
}
