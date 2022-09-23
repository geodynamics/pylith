#!/usr/bin/env nemesis
"""This is an extremely simple example showing how to set up an inversion
using PyLith-generated Green's functions. In this simple example,
there are no data uncertainties, and we use the minimum moment for a penalty function.

For the provided simulation in input files, simply run the script with no arguments.
The inversion results will be written to `output/step05_greensfns-inversion_results.txt`.

The order of the impulses and stations in the output depends on the parititioning among processes.
Consequently, we reorder the impulses and stations by the y coordinate to insure we maintain
consistent ordering.
"""

# Import argparse Python module (standard Python)
import argparse

# Import numpy and h5py modules (included in PyLith installation)
import numpy
import h5py

class InvertSlipApp:
    """Application to invert for fault slip using PyLith-generated Green's functions.
    """
    def __init__(self):
        self.filename_fault = "output/step05_greensfns-fault.h5"
        self.filename_responses = "output/step05_greensfns-gps_stations.h5"
        self.filename_observed = "output/step04_varslip-gps_stations.h5"
        self.filename_output = "output/step06_inversion-results.txt"
        self.penalties = [0.01, 0.1, 1.0]

    def main(self):
        """Entry point for running application.
        """
        args = self._parse_command_line()

        self.get_fault_impulses(args.filename_fault)
        self.get_station_responses(args.filename_responses)
        self.get_station_observed(args.filename_observed)

        self.penalties = list(map(float, args.penalties.split(",")))
        results = self.invert()
        self.write_results(self.filename_output, results)

    def get_fault_impulses(self, filename):
        """Get coordinates, amplitude, and order of impulses. Fault points are sorted by y coordinate.

        :param filename: Name of HDF5 file with fault data from PyLith Green's function simulation (Step 5).
        """
        h5 = h5py.File(filename, "r")
        y = h5['geometry/vertices'][:,1]
        slip = h5['vertex_fields/slip'][:,:,1]
        h5.close()

        # Sort fault points by y coordinate.
        reorder = numpy.argsort(y)

        self.impulse_y = y[reorder]
        self.impulse_slip = slip[:,reorder]

    def get_station_responses(self, filename):
        """
        Get coordinates and displacements at stations for Green's function responses.

        :param filename: Name of HDF5 file with point data from PyLith Green's function simulation (Step 5_).
        """
        h5 = h5py.File(filename, "r")
        xy = h5['geometry/vertices'][:,0:2]
        displacement = h5['vertex_fields/displacement'][:]
        h5.close()

        # Sort stations by y coordinate.
        reorder = numpy.argsort(xy[:,1])
        self.station_responses = displacement[:,reorder,:]

    def get_station_observed(self, filename):
        """
        Get coordinates and displacements at stations for fake observations.

        :param filename: Name of HDF5 file with point data from PyLith forward simulation (Step 4).
        """
        h5 = h5py.File(filename, "r")
        xy = h5['geometry/vertices'][:,0:2]
        displacement = h5['vertex_fields/displacement'][:]
        h5.close()

        # Sort stations by y coordinate.
        reorder = numpy.argsort(xy[:,1])
        self.station_observed = displacement[:,reorder,:]

    def invert(self):
        """Invert observations for amplitude of impulses and fault slip.
        """
        # Determine matrix sizes and set up A-matrix.
        nfaultpts = self.impulse_y.shape[0]
        nimpulses = self.station_responses.shape[0]
        nobs = self.station_responses.shape[1] * self.station_responses.shape[2]
        mat_A = self.station_responses.reshape((nimpulses, nobs)).transpose()

        # Data vector is a flattened version of the dataVals, plus the a priori
        # values of the parameters (assumed to be zero).
        vec_data = numpy.concatenate((self.station_observed.flatten(), 
            numpy.zeros(nimpulses, dtype=numpy.float64)))

        # Determine number of inversions and create empty array to hold results.
        ninversions = len(self.penalties)
        results = numpy.zeros((nfaultpts, 1 + ninversions))
        results[:,0] = self.impulse_y

        # Loop over number of inversions.
        for i, penalty in enumerate(self.penalties):

            # Scale diagonal by penalty parameter, and stack A-matrix with penalty matrix.
            mat_penalty = penalty * numpy.eye(nimpulses, dtype=numpy.float64)
            mat_design = numpy.vstack((mat_A, mat_penalty))

            # Form generalized inverse matrix.
            mat_gen_inverse = numpy.dot(numpy.linalg.inv(numpy.dot(mat_design.T, mat_design)), mat_design.T)

            # Solution is matrix-vector product of generalized inverse with data vector.
            impulse_amplitude = numpy.dot(mat_gen_inverse, vec_data)
            slip_soln = numpy.dot(impulse_amplitude, self.impulse_slip)
            results[:, 1+i] = slip_soln

            # Compute predicted results and residual.
            predicted = numpy.dot(mat_A, impulse_amplitude)
            residual = self.station_observed.flatten() - predicted
            residual_norm = numpy.linalg.norm(residual)
            print(f"Penalty parameter:  {penalty}")
            print(f"Residual norm:      {residual_norm}")
        return results

    def write_results(self, filename, results):
        """Write inversion results to file.

        :param filename: Name of output file.
        :param results: Inversion results as Numpy array.
        """
        header = "# y"
        for penalty in self.penalties:
            header += f" penalty={penalty}"
        header += "\n"

        with open(filename, "w") as fout:
            fout.write(header)
            numpy.savetxt(fout, results, fmt="%14.6e")

    def _parse_command_line(self):
        """Parse command line arguments.

        :returns: Command line arugment information as argparse.Namespace.
        """
        parser = argparse.ArgumentParser()

        parser.add_argument("--impulse-data", action="store", type=str, dest="filename_fault", default=self.filename_fault, help="Fault HDF5 data file from Green's function simulation.")
        parser.add_argument("--impulse-responses", action="store", type=str, dest="filename_responses", default=self.filename_responses, help="Station HDF5 data file from Green's function simulation.")
        parser.add_argument("--observed-data", action="store", type=str, dest="filename_observed", default=self.filename_observed, help="Station HDF5 data file from variable slip simulation.")
        parser.add_argument("--penalties", action="store", type=str, dest="penalties", default="0.01,0.1,1.0", help="Comma separated list of penalties.")
        parser.add_argument("--output", action="store", type=str, dest="filename_output", default=self.filename_output, help="Name of output file with inversion results.")
        return parser.parse_args()

if __name__ == "__main__":
    InvertSlipApp().main()


# End of file
