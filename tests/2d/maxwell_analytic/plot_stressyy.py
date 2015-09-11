#!/usr/bin/env python


def plot(sim):
    import h5py
    import pylab
    import numpy

    filename = "output/%s-statevars.h5" % sim
    h5 = h5py.File(filename, "r")
    stress = h5['cell_fields/stress'][:,:,1]
    t = h5['time'][:].ravel()
    h5.close()

    filename = "output/%s-statevars_info.h5" % sim
    h5 = h5py.File(filename, "r")
    mat_mu = h5['cell_fields/mu'][0,0,0]
    mat_lambda = h5['cell_fields/lambda'][0,0,0]
    mat_density = h5['cell_fields/density'][0,0,0]
    mat_tm = h5['cell_fields/maxwell_time'][0,0,0]
    h5.close()


    K = mat_lambda + 2.0/3.0*mat_mu
    G = mat_mu
    viscosity = mat_tm * mat_mu
    theta = viscosity / G
    analytic = -10e+6*(1.0-6*G/(3*K+4*G)*numpy.exp(-3*K*t/((3*K+4*G)*theta)))
    
    pylab.plot(t, analytic[:], 'k-', t, stress[:,0], 'r--')
    pylab.show()

    return

# ======================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--sim", action="store", dest="sim", required=True)
    args = parser.parse_args()

    plot(args.sim)


