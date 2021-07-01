#!/usr/bin/env nemesis


def plot(sim):
    import h5py
    import pylab
    import numpy

    filename = "output/%s-maxwell.h5" % sim
    h5 = h5py.File(filename, "r")
    stress = h5['vertex_fields/cauchy_stress'][:, :, 1]
    t = h5['time'][:].ravel()
    h5.close()

    filename = "output/%s-maxwell_info.h5" % sim
    h5 = h5py.File(filename, "r")
    shear_modulus = h5['vertex_fields/shear_modulus'][0, 0, 0]
    bulk_modulus = h5['vertex_fields/bulk_modulus'][0, 0, 0]
    maxwell_time = h5['vertex_fields/maxwell_time'][0, 0, 0]
    h5.close()

    viscosity = maxwell_time * shear_modulus
    theta = viscosity / shear_modulus
    analytic = -10e+6*(1.0-6*shear_modulus/(3*bulk_modulus+4*shear_modulus)
                       * numpy.exp(-3*bulk_modulus*t/((3*bulk_modulus+4*shear_modulus)*theta)))

    pylab.plot(t, analytic[:], 'k-', t, stress[:, 0], 'r--')
    pylab.show()

    return


# ======================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--sim", action="store", dest="sim", required=True)
    args = parser.parse_args()

    plot(args.sim)
