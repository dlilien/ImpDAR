/*
 * mig_cython.c
 * Copyright (C) 2019 dlilien <dlilien@berens>
 *
 * Distributed under terms of the MIT license.
 */

#include "mig_cython.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <time.h> 
#include <signal.h>
#include <stdlib.h>
#include "interpolation_c.h"

void mig1k (float k, float fmax, float speed, int nt,
        float dt, float *v, double complex *p, double complex *q);
static void closefiles(void);

/*  Migrate data into migdata */
void mig_kirch_loop (double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, int nearfield){
    int i, j, k, l;
    double min, m;
    double costheta;
    double integral;
    double rs;
    int Didx;
    int lmin_prev, rmin_prev;
    clock_t iter_start, this_time;
    float time_remaining, el_t;
    setbuf(stdout, NULL);
    printf("Iter: ");

    
    iter_start = clock();
    /* You can flip the loop order, but I think this is marginally faster
     * based on c row-major ordering. Could be wrong though... */
    for(j=0;j<tnum;j++){
        if (j % 100 == 0){
            printf("%d", j);
        }else{
            printf(".");
        }
        for(i=0;i<snum;i++){
            integral = 0.0;
            rmin_prev = i;
            lmin_prev = i;

            /* Do left and right separately so we know when to break
            * the break lets us skip a lot of iterations for long profiles */
            for(k=j + 1;k<tnum;k++){
                /* rs is the distance to the closest point in this trace??? */
                rs = sqrt(pow(dist[k] - dist[j], 2.) + zs2[i]);
                /* We do not get closer as we go, so we can break the loop
                 * if we are far away */
                if(2. * rs / vel > max_travel_time){
                     break;
                }

                costheta = zs[i] / rs;
                min = 1.0e6;
                /* The new min should not be smaller than the previous rmin
                * since otherwise the hyperbola would be inverted (smiles) */
                for(l=rmin_prev;l<snum;l++){
                    m = fabs(tt_sec[l] - 2. * rs / vel);
                    if(m < min){
                        min = m;
                        Didx = l;
                    }
                }
                rmin_prev = Didx;
                integral += gradD[Didx * tnum + k] * costheta / vel;
            }

            for(k=j;k>=0;k--){
                rs = sqrt(pow(dist[k] - dist[j], 2.) + zs2[i]);
                if(2. * rs / vel > max_travel_time){
                     break;
                }
                costheta = zs[i] / rs;
                if (rs <= 1.0e8){
                    costheta = 0.;
                }
                min = 1.0e6;
                for(l=lmin_prev;l<snum;l++){
                    m = fabs(tt_sec[l] - 2. * rs / vel);
                    if(m < min){
                        min = m;
                        Didx = l;
                    }
                }
                lmin_prev = Didx;
                integral += gradD[Didx * tnum + k] * costheta / vel;
            }

            migdata[i * tnum + j] = integral;
        }
        if (j % 100 == 0){
            if (j > 0){
                this_time = clock();
                el_t = (this_time - iter_start) / CLOCKS_PER_SEC;
                time_remaining = (tnum - j) * el_t / j;
                printf("\nEst. time remaining: %f sec", time_remaining);
            }
        }
    }
    printf("\n");
}


/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */
/* SUMIGTK: $Revision: 1.17 $ ; $Date: 2011/11/16 22:14:43 $		*/
/* Modification 08/2019 by D. Lilien to make python callable */
/* Credits:
 *	CWP: Dave Hale
 */

char *sdoc[] = {
"									",
" SUMIGTK - MIGration via T-K domain method for common-midpoint stacked data",
"									",
" sumigtk <stdin >stdout dxcdp= [optional parms]			",
"									",
" Required Parameters:							",
" dxcdp                   distance between successive cdps		",
"									",
" Optional Parameters:							",
" fmax=Nyquist            maximum frequency				",
" tmig=0.0                times corresponding to interval velocities in vmig",
" vmig=1500.0             interval velocities corresponding to times in tmig",
" vfile=                  binary (non-ascii) file containing velocities v(t)",
" nxpad=0                 number of cdps to pad with zeros before FFT	",
" ltaper=0                length of linear taper for left and right edges", 
" verbose=0               =1 for diagnostic print			",
"									",
" Notes:								",
" Input traces must be sorted by either increasing or decreasing cdp.	",
" 									",
" The tmig and vmig arrays specify an interval velocity function of time.",
" Linear interpolation and constant extrapolation is used to determine	",
" interval velocities at times not specified.  Values specified in tmig	",
" must increase monotonically.						",
"									",
" Alternatively, interval velocities may be stored in a binary file	",
" containing one velocity for every time sample.  If vfile is specified,",
" then the tmig and vmig arrays are ignored.				",
"									",
" The time of first sample is assumed to be zero, regardless of the value",
" of the trace header field delrt.					",
" 									",
" The migration is a reverse time migration in the (t,k) domain. In the	",
" first step, the data g(t,x) are Fourier transformed x->k into the	",	
" the time-wavenumber domain g(t,k).					",
"									",
" Then looping over wavenumbers, the data are then reverse-time		",
" finite-difference migrated, wavenumber by wavenumber.  The resulting	",
" migrated data m(tau,k), now in the tau (migrated time) and k domain,	",
" are inverse fourier transformed back into m(tau,xout) and written out.",	
"									",
NULL};


void sumigtk(double * migdata, double * data, int snum, int tnum, float * vt, double speed, double dxcdp, double fmax, int tnumpad, int ltaper, int verbose){
	int tnumfft,ix,it,nk,ik,
		ntmig,nvmig,itmig,i2;
	float dt,dx,dk,taper,t,k,fftscl,
		*tmig,*vmig,**gtx;
	double complex **gtk;
	char *vfile="";

	dt = dt / 1000000.0;

	if (verbose) fprintf(stderr,"\t%d traces input\n",tnum);
	
	/* determine wavenumber sampling */
	tnumfft = npfaro(tnum+tnumpad,2*(tnum+tnumpad));
	nk = tnumfft/2+1;
	dk = 2.0*PI/(tnumfft*dx);

	/* allocate space for Fourier transform */
    gtk = (double complex**)malloc(snum*sizeof(void*));
	for (i2=0; i2<nk; i2++){
		gtk[i2] = (double complex*)gtk[0]+sizeof(double complex)*snum*i2;
    }

	gtx = (float**)(tnumfft, sizeof(void*));
	for (ix=0; ix<tnumfft; ++ix)
		gtx[ix] = (float*)gtk[0]+ix*snum;

	/* read and apply fft scaling to traces and pad with zeros */
	fftscl = 1.0/tnumfft;
	for (ix=0; ix<tnum; ++ix) {
		efread(gtx[ix],sizeof(float),snum,tracefp);
		for (it=0; it<snum; ++it)
			gtx[ix][it] *= fftscl;
		if (ix<ltaper) {
			taper = (float)(ix+1)/(float)(ltaper+1);
			for (it=0; it<snum; ++it)
				gtx[ix][it] *= taper;
		} else if (ix>=tnum-ltaper) {
			taper = (float)(tnum-ix)/(float)(ltaper+1);
			for (it=0; it<snum; ++it)
				gtx[ix][it] *= taper;
		}
	}
	for (ix=tnum; ix<tnumfft; ++ix)
		for (it=0; it<snum; ++it)
			gtx[ix][it] = 0.0;
	
	/* Fourier transform g(t,x) to g(t,k) */
	pfa2rc(-1,2,snum,tnumfft,gtx[0],gtk[0]);
	if (verbose) fprintf(stderr,"\tFourier transform done\n");
	
	/* loop over wavenumbers */
	for (ik=0,k=0.0; ik<nk; ++ik,k+=dk) {
	
		/* report */
		if (verbose && ik%(nk/10>0?nk/10:1)==0)
			fprintf(stderr,"\t%d of %d wavenumbers done\n",
				ik,nk);
		
		/* migrate */
		mig1k(k,fmax,speed,snum,dt,vt,gtk[ik],gtk[ik]);
	}
	
	/* Fourier transform g(t,k) to g(t,x) */
	pfa2cr(1,2,snum,tnumfft,gtk[0],gtx[0]);
	if (verbose) fprintf(stderr,"\tinverse Fourier transform done\n");
	
	/* output migrated traces with headers */
	for (ix=0; ix<tnum; ++ix) {
		efread(&tr,HDRBYTES,1,headerfp);
		memcpy( (void *) tr.data,
				(const void *) gtx[ix], snum*sizeof(float));
		puttr(&tr);
	}
}

void mig1k (float k, float fmax, float speed, int nt, float dt, 
	float *v, double complex *p, double complex *q)
/*****************************************************************************
migration in t-k domain for one wavenumber
******************************************************************************
Input:
k		wavenumber
fmax		maximum frequency
speed		speed parameter - >>1.0 for lots of dispersion
nt		number of time samples
dt		time sampling interval
v		array[nt] containing interval velocities v(t)
p		array[nt] containing input p(t;k)

Output:
q		array[nt] containing migrated q(t;k)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/05/90
*****************************************************************************/
{
	int nfac=5,ns,is,it,i1,i2,i1stop;
	float fac=0.75,ds,g00r,g00i,g01r,g01i,g10r,g10i,g11r,g11i,
		temp,vmin,vmax,ct,
		*vs,*cs,*cp,*grp,*gip;
	double complex czero=cmplx(0.0,0.0),*gs;
	
	/* determine time sampling to avoid excessive grid dispersion */
	for (it=1,vmin=vmax=v[0]; it<nt; ++it) {
		if (v[it]<vmin) vmin = v[it];
		if (v[it]>vmax) vmax = v[it];
	}
	if (k!=0.0) {
		ds = fac*MAX(0.4*PI/ABS(vmax*k),0.1*vmin/(vmax*fmax));
		ds *= speed;
		ds = MIN(ds,dt);
	} else {
		ds = dt;
	}
	ns = 1+(nt-1)*dt/ds;
	ns = MIN(ns,nfac*nt);
	ds = (nt-1)*dt/(ns-1);
	fprintf(stderr,"ns=%d\n",ns);
	
	/* allocate workspace */
	vs = (float *)calloc(ns, sizeof(float));
	cs = (float *)calloc(ns*2, sizeof(float));
	gs = (double complex *)calloc(ns, sizeof(double complex));
	
	/* resample v(t) and p(t) */
	ress8r(nt,dt,0.0,v,v[0],v[nt-1],ns,ds,0.0,vs);
	ress8c(nt,dt,0.0,p,czero,czero,ns,ds,0.0,gs);

        /* compute finite-difference coefficients */
	for (is=0; is<ns; is++) {
		temp = 0.125*vs[is]*k*ds;
		temp = temp*temp;
		temp = (1.0-temp)/(1.0+temp);
		cs[2*is] = cs[2*is+1] = temp;
	}

	/* loop over t2 = (tau-t)/sqrt(2) */
	for (i2=2-ns; i2<ns; i2++) {

		/* determine t1 stop index */
		i1stop = (i2<=0)?1-i2:i2;

		/* initialize finite-difference star and coefficient */
		g00r = g00i = g01r = g01i = 0.0;
		grp = (float*)(&gs[ns-1]);
		gip = grp+1;
		cp = &cs[i2+ns-2];
		ct = *cp--;

		/* loop over t1 = (tau+t)/sqrt(2) */
		for (i1=ns-1; i1>=i1stop; i1--) {

			/* update real part */
			g10r = g00r;
			g11r = g01r;
			g00r = *grp;
			*grp = g01r = ct*(g11r+g00r)-g10r;

			/* update imaginary part */
			g10i = g00i;
			g11i = g01i;
			g00i = *gip;
			*gip = g01i = ct*(g11i+g00i)-g10i;

			/* update pointers and finite-difference coefficient */
			grp -= 2; gip -= 2;
			ct = *cp--;
		}
	}
	
	/* resample q(t) */
	ress8c(ns,ds,0.0,gs,czero,czero,nt,dt,0.0,q);
	
	/* free workspace */
	free(vs);
	free(cs);
	free(gs);
}
