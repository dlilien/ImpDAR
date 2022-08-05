/*
 * coherence.c
 * Copyright (C) 2021 David Lilien <david.lilien@umanitoba.ca>
 *
 * Distributed under terms of the GNU GPL3.0 license.
 */

#include <complex.h>
#include "coherence.h"
#include <math.h>
#include <stdio.h>
#include <time.h> 
#include <signal.h>
#include <stdlib.h>
#include <Python.h>

void coherence2d (double complex * chhvv, double complex * HH, double complex * VV, int nrange, int ntheta, int range_bins, int azimuth_bins){
    int i, j, ii, jj;
    int imin, imax, jmin, jmax;
    double complex numerator;
    double complex denominator1;
    double complex denominator2;
    PySys_WriteStdout("Beginning iteration through %d azimuths...\nAzimuth bin: ", azimuth_bins);
    fflush(stdout);

    /* You can flip the loop order, but I think this is marginally faster
     * based on c row-major ordering. Could be wrong though... */

    for(i=ntheta;i<azimuth_bins - ntheta;i++){
        if (i % 10 == 0){
            if (i > 0){
                PySys_WriteStdout("\n%d", i);
                fflush(stdout);
            }
        } else {
                PySys_WriteStdout(".");
                fflush(stdout);
        }

        for(j=0;j<range_bins;j++){
            imin = i - ntheta;
            imax = i + ntheta;
            jmin = max(0, j - nrange);
            jmax = min(range_bins - 1, j + nrange);
            numerator = 0.0 + 0.0 * _Complex_I;
            denominator1 = 0.0 + 0.0 * _Complex_I;
            denominator2 = 0.0 + 0.0 * _Complex_I;
            for (ii=imin; ii<imax; ii++){
                for (jj=jmin; jj<jmax; jj++){
                    numerator += HH[jj * azimuth_bins + ii] * conj(VV[jj * azimuth_bins + ii]);
                    denominator1 += pow(cabs(HH[jj * azimuth_bins + ii]), 2); 
                    denominator2 += pow(cabs(VV[jj * azimuth_bins + ii]), 2);
                }
            }
            chhvv[j * azimuth_bins + i] = numerator/(sqrt(denominator1) * sqrt(denominator2));
        }
    }
}
