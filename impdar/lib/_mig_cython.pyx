""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """
import sys
cimport numpy as np
np.import_array()
import numpy as np

# cdefine the signature of our c function
# cdef extern from "mig_cython.h":
#    void mig_cython (double * data, double * migdata, int tnum, int snum, double * dist, double * zs, double * zs2, double * tt_sec, double vel, double * gradD, double max_travel_time, bint nearfield)

# create the wrapper code, with numpy type annotations
def migrationKirchoffLoop(np.ndarray[double, ndim=2, mode="c"] data not None,
                          np.ndarray[double, ndim=2, mode="c"] migdata not None,
                          int tnum,
                          int snum,
                          np.ndarray[double, ndim=1, mode="c"] dist not None,
                          np.ndarray[double, ndim=1, mode="c"] zs not None,
                          np.ndarray[double, ndim=1, mode="c"] zs2 not None,
                          np.ndarray[double, ndim=1, mode="c"] tt_sec not None,
                          float vel,
                          np.ndarray[double, ndim=2, mode="c"] gradD not None,
                          float max_travel_time,
                          bint nearfield
                          ):
    """
    mig_cython(<double*> np.PyArray_DATA(data),
               <double*> np.PyArray_DATA(migdata),
               tnum,
               snum,
               <double*> np.PyArray_DATA(dist),
               <double*> np.PyArray_DATA(zs),
               <double*> np.PyArray_DATA(zs2),
               <double*> np.PyArray_DATA(tt_sec),
               vel,
               <double*> np.PyArray_DATA(gradD),
               max_travel_time,
               nearfield
               )
               """
            # Loop through all traces

    cdef int xi
    cdef int ti
    cdef float x
    cdef np.ndarray dists2
    cdef np.ndarray rs
    cdef np.ndarray gradDhyp
    cdef np.ndarray Dhyp
    cdef float integral
    cdef np.ndarray Didx

    print('Migrating trace number:')
    for xi in range(tnum):
        print('{:d}, '.format(xi))
        sys.stdout.flush()
        # get the trace distance
        x = dist[xi]
        dists2 = (dist - x)**2.
        # Loop through all samples
        for ti in range(snum):
            # get the radial distances between input point and output point
            rs = np.sqrt(dists2 + zs2[ti])
            # find the cosine of the angle of the tangent line, correct for obliquity factor
            with np.errstate(invalid='ignore'):
                costheta = zs[ti] / rs
            # get the exact indices from the array (closest to rs)
            Didx = np.argmin(np.abs(np.atleast_2d(tt_sec).transpose() - 2. * rs / vel), axis=0)
            # integrate the farfield term
            gradDhyp = gradD[Didx, np.arange(len(Didx))]
            gradDhyp[2. * rs / vel > max_travel_time] = 0.    # zero points that are outside of the domain
            integral = np.nansum(gradDhyp * costheta / vel)  # TODO: Yilmaz eqn 4.5 has an extra r in this weight factor???
            # integrate the nearfield term
            if nearfield:
                Dhyp = data[Didx, np.arange(len(Didx))]
                Dhyp[2. * rs / vel > max_travel_time] = 0.    # zero points that are outside of the domain
                integral += np.nansum(Dhyp * costheta / rs**2.)
            # sum the integrals and output
            migdata[ti, xi] = 1. / (2. * np.pi) * integral
