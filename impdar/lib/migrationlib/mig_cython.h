/*
 * mig_cython.h
 * Copyright (C) 2019 dlilien <dlilien@berens>
 *
 * Distributed under terms of the GNU GPL3.0 license.
 */

#ifndef MIG_CYTHON_H
#define MIG_CYTHON_H

void mig_kirch_loop (double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, int nearfield);
void sumigtk(double * migdata, double * data, int ns, int nt, float * vt, double speed, double dxcdp, double fmax, int nxpad, int ltaper, int verbose);

#endif /* !MIG_CYTHON_H */
