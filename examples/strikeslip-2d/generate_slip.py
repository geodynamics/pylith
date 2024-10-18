#!/usr/bin/env nemesis
"""Python script for generating slip profile.

REQUIRES: scipy (not included in PyLith binary distribution)
"""

import math
import abc

import scipy
import numpy

from pythia.pyre.units.length import m, km


def cli():
    """Command line interface.
    """
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--eq-magnitude", action="store", type=float, dest="eq_magnitude", default=6.5, help="Earthquake magnitude.")
    parser.add_argument("--total-length", action="store", type=float, dest="total_length", default=40.0, help="Rupture length (km).")
    parser.add_argument("--taper-length", action="store", type=float, dest="taper_length", default=0.5, help="Taper length (km).")
    parser.add_argument("--resolution", action="store", type=float, dest="resolution", default=100.0, help="Resolution (m).")
    parser.add_argument("--db-filename", action="store", dest="filename_db", default="slip_variable.spatialdb", help="Filename for slip spatial database.")
    parser.add_argument("--raw-filename", action="store", dest="filename_raw", default="output/slip_variable.txt", help="Filename for raw slip data.")
    args = parser.parse_args()

    generator = SlipGenerator(
        eq_magnitude=args.eq_magnitude,
        total_length=args.total_length*km,
        taper_length=args.taper_length*km,
        dx=args.resolution*m,
        )
    (x, slip) = generator.generate()
    SpatialDBWriter.write(args.filename_db, x, slip)
    RawWriter.write(args.filename_raw, x, slip)


class SlipGenerator:
    """Object for generating 1D slip profile.
    """

    def __init__(self, eq_magnitude: float, total_length: float, taper_length: float, dx: float):
        """Generate 1D slip profile.

        Args:
            eq_magnitude: Earthquake magnitude.
            total_length: Rupture length (m)
            taper_length: Length of taper (m) at ends of rupture.
            dx: Discretization size for slip profile.
        """
        self.eq_magnitude = eq_magnitude
        self.total_length = total_length
        self.taper_length = taper_length
        self.dx = dx

    def generate(self) -> tuple:
        """Generate slip profile.

        Returns:
            (numpy.ndarray, numpy.ndarray): (distance along fault (m), slip along fault (m))
        """
        MAX_ITERATIONS = 100
        END_TOLERANCE = 0.02
        ZERO_TOLERANCE = 0.05
        SPLINE_SPACING = 4 # Ratio of discretization size in raw slip profile to discretization size in spline

        dx_points = SPLINE_SPACING*self.dx # Discretization size of raw slip profile
        randomizer = RandomizerLA(
            eq_magnitude=self.eq_magnitude, 
            total_length=self.total_length,
            dx=dx_points,
        )
        total_length_m = self.total_length.value
        dx_points_m = dx_points.value
        x_pts = numpy.arange(-0.5*total_length_m, +0.5*(total_length_m+dx_points_m), dx_points_m)
        for i in range(MAX_ITERATIONS):
            # Iterate until we find slip profile that is nearly zero at the ends
            slip_pts = randomizer.create(x_pts)
            if slip_pts[0] < END_TOLERANCE and numpy.sum(slip_pts == 0.0) < (1.0-ZERO_TOLERANCE)*len(slip_pts):
                break
        else:
            print("Failed to generated slip distribution.")

        # Apply taper to force slip profile to zero at the ends.
        taper = Taper(
            total_length=self.total_length,
            taper_length=self.taper_length,
        )
        taper.apply(x_pts, slip_pts)

        # Fit spline to slip profile
        dx_spline_m = self.dx.value
        x = numpy.arange(-0.5*total_length_m, +0.5*(total_length_m+dx_spline_m), dx_spline_m)
        spline = scipy.interpolate.CubicSpline(x_pts, slip_pts)
        slip_spline = spline(x)
        slip_spline[0] = 0.0
        slip_spline[-1] = 0.0
        slip_spline[slip_spline < 0.0] = 0.0
        return (x, slip_spline)


class SpatialDBWriter:
    """Object for writing rupture information to spatial database file.
    """

    @staticmethod
    def write(filename: str, y: numpy.ndarray, slip: numpy.ndarray):
        """Write rupture information to SimpleGridDB ASCII file.

        Args:
            filename: Name of spatial database file.
            y: y coordinate (m)
            slip: Slip (m) along fault
        """
        from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
        from spatialdata.geocoords.CSCart import CSCart

        DOMAIN_Y = 150.0*km

        writer = SimpleGridAscii()
        writer.filename = filename

        # Add points to encompass entire fault (regions beyond rupture)
        y_ext = numpy.concatenate(([-0.5*DOMAIN_Y.value], y, [+0.5*DOMAIN_Y.value]))
        slip_ext = numpy.concatenate(([0.0], slip, [0.0]))

        cs = CSCart()
        cs.spaceDim = 2
        cs._configure()
        points = numpy.zeros((len(y_ext), 2))
        points[:,1] = y_ext

        data = {
            "points": points,
            'x': numpy.array([0]),
            'y': y_ext,
            'z': numpy.array([0]),
            "coordsys": cs,
            "data_dim": 1,
            "values": [{
                "name": "final_slip_left_lateral",
                "units": "m",
                "data": slip_ext,
            },{
                "name": "final_slip_opening",
                "units": "m",
                "data": 0*slip_ext,
            },{
                "name": "initiation_time",
                "units": "s",
                "data": 0*slip_ext,
            }]
        }
        writer.write(data)


class RawWriter:
    """Write slip profile to ASCII file for easy plotting.
    """

    @staticmethod
    def write(filename: str, x: numpy.ndarray, slip: numpy.ndarray):
        """Write rupture information to SimpleGridDB ASCII file.

        Args:
            filename: Name of spatial database file.
            x: Along-fault coordinate (m)
            slip: Slip (m) along fault
        """
        data = numpy.stack((x, slip)).T
        numpy.savetxt(filename, data, fmt=("%8.1f", "%5.3f"))


class Randomizer(abc.ABC):
    """Abstract base class for generating random distribution.
    """
    SHEAR_MODULUS = 3.3e+10 # Pa
    SEED = 234

    def __init__(self, eq_magnitude: float, total_length: float, dx: float):
        """Constructor.

        Args:
            eq_magnitude: Earthquake magnitude.
            total_length: Total length of rupture.
            dx: Discretization size for slip profile along fault.
        """
        self.eq_magnitude = eq_magnitude
        self.total_length = total_length
        self.dx = dx
        numpy.random.seed(self.SEED)

    @staticmethod
    def _inverse_fft(amplitude: numpy.ndarray, phase: numpy.ndarray) -> numpy.ndarray:
        """Compute inverse FFT to get slip profile from amplitude and phase.

        Args:
            amplitude: Slip amplitude.
            phase: Slip phase.

        Returns:
            Slip profile.
        """
        amplitudeFFT = amplitude * (numpy.cos(phase) + numpy.sin(phase)*complex(0,1))
        return numpy.real(numpy.fft.ifft(amplitudeFFT))

    @staticmethod
    def _avg_slip(eq_magnitude: float) -> float:
        """Compute average slip given earthquake magnitude.

        Args:
            eq_magnitude: Earthquake magnitude.

        Returns:
            Average slip (m).
        """
        area = 10**(0.8*(eq_magnitude - 3.30))
        area *= 1.0e+6  # km**2 -> m**2
        seismic_moment = 10**(1.5*(eq_magnitude+10.7)-7)
        avg_slip = seismic_moment / (area * Randomizer.SHEAR_MODULUS)
        return avg_slip
        
    @staticmethod
    def _correlation_coef(eq_magnitude: float) -> float:
        """Compute correlation coefficient.

        Args:
            eq_magnitude: Earthquake magnitude.

        Returns:
            Correlation coefficient (m).
        """
        return 1*km * 10**(0.5*eq_magnitude - 2.5)


class RandomizerLA(Randomizer):
    """Random slip profile based on Lavrentiadis and Abrahamson
    (2019, doi: 10.1785/0120180252).

    We use a slightly larger value for phi_0.
    """
    
    Np = 1.24
    #PHI0 = 0.265 # Original value used by Lavrentiadis and Abrahamson
    PHI0 = 0.35 # Modified value

    def create(self, x: numpy.ndarray) -> numpy.ndarray:
        """Create 1D slip profile.

        Args:
            x: Along-fault coordinate.

        Returns:
            Slip profile.
        """
        mag = self._create_amplitude(x)
        phase = self._create_phase(x)

        slip = numpy.abs(self._inverse_fft(mag, phase))
        i_min = numpy.argmin(slip)
        slip = numpy.concatenate((slip[i_min:], slip[:i_min]))
        return slip

    def _create_amplitude(self, x: numpy.ndarray) -> numpy.ndarray:
        """Generate slip amplitude.

        Args:
            x: Along-fault coordinate.

        Returns:
            Slip amplitude.
        """
        avg_slip_m = self._avg_slip(self.eq_magnitude)

        kx = numpy.fft.fftfreq(len(x), self.dx / km)
        k = numpy.abs(kx)
        kc = self._kc()
        scale0 = avg_slip_m * len(x)
        scale = scale0 / numpy.sqrt(1.0+(k/kc)**(2*self.Np))

        amplitude = scale * 10**(numpy.random.normal(loc=0.0, scale=self.PHI0, size=len(scale)))
        amplitude[k==0] = scale0
        return amplitude

    def _create_phase(self, x: numpy.ndarray) -> numpy.ndarray:
        """Generate slip phase.

        Args:
            x: Along-fault coordinate.

        Returns:
            Slip phase.
        """
        total_length_km = self.total_length / km
        dx_km = self.dx / km

        mu = -total_length_km * math.pi
        phase_der = numpy.random.logistic(loc=mu, scale=self._s(), size=len(x))

        dk = 1.0 / (dx_km * len(x))
        phase = scipy.integrate.cumulative_simpson(phase_der, dx=dk, initial=0.0)
        phase = phase % (2.0*math.pi)
        return phase

    def _kc(self) -> float:
        """Compute kc parameter.

        Returns:
            kc parameter.
        """
        c1 = -2.03
        c2 = -1.0
        c3 = 1.6
        total_length_km = self.total_length / km
        return 10**(c1 + c2*(math.log10(total_length_km) - c3))

    def _s(self) -> float:
        """Compute asa parameter.

        Returns:
            `s` parameter.
        """
        d1 = 1.49
        d2 = 1.0
        d3 = 1.6
        total_length_km = self.total_length / km
        return 10**(d1 + d2*(math.log10(total_length_km) - d3))


class Taper:
    """Object for applying linear taper at ends of slip profile.
    """

    def __init__(self, total_length: float, taper_length: float):
        """Constructor.

        Args:
            total_length: Total rupture length.
            taper_length: Length of linear taper at each end.
        """
        self.total_length = total_length
        self.taper_length = taper_length

    def apply(self, x: numpy.ndarray, slip: numpy.ndarray):
        """Apply linear taper at ends of slip profile.

        Args:
            x: Along-fault coordinate.
            slip: Slip profile to apply taper to.
        """
        l_taper = self.taper_length.value
        maskL = x <= x[0] + l_taper
        maskR = x >= x[-1] - l_taper
        taper = numpy.logical_and(~maskL, ~maskR) + maskL * (x-x[0]) / l_taper + maskR * (x[-1] - x) / l_taper
        slip *= taper


if __name__ == "__main__":
    cli()
