/*
 * coherence.h
 * Copyright (C) 2021 David Lilien <david.lilien@umanitoba.ca>
 *
 * Distributed under terms of the GNU GPL3.0 license.
 */
#include <complex.h>

#ifndef COHERENCE_h
#define COHERENCE_H
#define min(X,Y) (((X) < (Y)) ? (X) : (Y))
#define max(X,Y) (((X) > (Y)) ? (X) : (Y))
void coherence2d (double complex * chhvv, double complex * HH, double complex * VV, int nrange, int ntheta, int range_bins, int azimuth_bins);
#endif /* !COHERENCE_H */
