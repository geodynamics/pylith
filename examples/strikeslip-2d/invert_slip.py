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

import pathlib

# Import numpy and h5py modules (included in PyLith installation)
import numpy
import h5py

OUTPUT_DIR = pathlib.Path("output")

def cli():
    """Command line interface.
    """
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--impulse-data", action="store", type=str, dest="filename_fault", default="step05_greensfns-fault.h5", help="Fault HDF5 data file from Green's function simulation.")
    parser.add_argument("--impulse-responses", action="store", type=str, dest="filename_responses", default="step05_greensfns-gnss_stations.h5", help="Station HDF5 data file from Green's function simulation.")
    parser.add_argument("--observed-data", action="store", type=str, dest="filename_observed", default="step04_varslip-gnss_stations.h5", help="Station HDF5 data file from variable slip simulation.")
    parser.add_argument("--penalties", action="store", type=str, dest="penalties", default="0.02,0.1,0.4", help="Comma separated list of penalties.")
    parser.add_argument("--output", action="store", type=str, dest="filename_output", default="step06_inversion-results.txt", help="Name of output file with inversion results.")

    args = parser.parse_args()
    app = InvertSlipApp(output_dir=OUTPUT_DIR)
    app.run(
        filename_fault=args.filename_fault,
        filename_responses=args.filename_responses,
        filename_observed=args.filename_observed,
        filename_output=args.filename_output,
        penalties=[float(p) for p in args.penalties.split(",")],
        )


class InvertSlipApp:
    """Application to invert for fault slip using PyLith-generated Green's functions.
    """
    def __init__(self, output_dir: str):
        """Constructor

        Args:
            output_dir: Directory with simulation output.
        """
        self.output_dir = output_dir

    def run(self, filename_fault: str, filename_responses: str, filename_observed: str, filename_output: str, penalties: list):
        """Entry point for running application.
        """
        self._get_fault_impulses(self.output_dir / filename_fault)
        self._get_station_responses(self.output_dir / filename_responses)
        self._get_station_observed(self.output_dir / filename_observed)

        self.penalties = penalties
        results = self._invert()
        self.write_results(self.output_dir / filename_output, results)

    def _get_fault_impulses(self, filename: str):
        """Get coordinates, amplitude, and order of impulses. Fault points are sorted by y coordinate.

        Args:
            filename: Name of HDF5 file with fault data from PyLith Green's function simulation (Step 5).
        """
        h5 = h5py.File(filename, "r")
        y = h5['geometry/vertices'][:,1]
        slip = h5['vertex_fields/slip'][:,:,1]
        h5.close()

        # Sort fault points by y coordinate.
        reorder = numpy.argsort(y)

        self.impulse_y = y[reorder]
        self.impulse_slip = slip[:,reorder]

    def _get_station_responses(self, filename):
        """
        Get coordinates and displacements at stations for Green's function responses.

        Args:
            filename: Name of HDF5 file with point data from PyLith Green's function simulation (Step 5_).
        """
        h5 = h5py.File(filename, "r")
        xy = h5['geometry/vertices'][:,0:2]
        displacement = h5['vertex_fields/displacement'][:]
        h5.close()

        # Sort stations by y coordinate.
        reorder = numpy.argsort(xy[:,1])
        self.station_responses = displacement[:,reorder,:]

    def _get_station_observed(self, filename):
        """
        Get coordinates and displacements at stations for fake observations.

        Args:
            filename: Name of HDF5 file with point data from PyLith forward simulation (Step 4).
        """
        h5 = h5py.File(filename, "r")
        xy = h5['geometry/vertices'][:,0:2]
        displacement = h5['vertex_fields/displacement'][:]
        h5.close()

        # Sort stations by y coordinate.
        reorder = numpy.argsort(xy[:,1])
        self.station_observed = displacement[:,reorder,:]

    def _invert(self) -> numpy.ndarray:
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

    def write_results(self, filename: str, results: numpy.ndarray):
        """Write inversion results to file.

        Args:
            filename: Name of output file.
            results: Inversion results as Numpy array.
        """
        header = "# y"
        for penalty in self.penalties:
            header += f" penalty={penalty}"
        header += "\n"

        with open(filename, "w") as fout:
            fout.write(header)
            numpy.savetxt(fout, results, fmt="%14.6e")


if __name__ == "__main__":
    cli()


# End of file
