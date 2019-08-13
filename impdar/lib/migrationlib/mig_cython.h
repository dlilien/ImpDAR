/*
 * mig_cython.h
 * Copyright (C) 2019 dlilien <dlilien@berens>
 *
 * Distributed under terms of the GNU GPL3.0 license.
 */

#ifndef MIG_CYTHON_H
#define MIG_CYTHON_H
#include <stdbool.h> 

void mig_kirch_loop (double * data, double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, bool nearfield, double * dist2, double * rs);
#endif /* !MIG_CYTHON_H */
