#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

# @file euler/poleconvert

# @brief Python application to convert Euler poles from geographic to
# Cartesian or vice-versa.

import math
import numpy

from pythia.pyre.applications.Script import Script as Application


class PoleConvert(Application):
    """Python application to convert Euler poles from geographic to Cartesian.
    """

    class Inventory(Application.Inventory):
        """Python object for managing PoleConvert facilities and properties.
        """

        # @class Inventory
        # Python object for managing PoleConvert facilities and properties.
        ##
        # \b Properties
        # @li \b poles_input_file Filename of file containing input poles.
        # @li \b poles_output_file Filename of output poles file.
        # @li \b geog_to_cart Flag to indicate geographic to Cartesian.
        ##

        import pythia.pyre.inventory

        polesInputFile = pythia.pyre.inventory.str("poles_input_file",
                                                   default="input.poles")
        polesInputFile.meta['tip'] = "Filename of file containing input poles."

        polesOutputFile = pythia.pyre.inventory.str("poles_output_file",
                                                    default="output.poles")
        polesOutputFile.meta['tip'] = "Filename of output poles file."

        geogToCart = pythia.pyre.inventory.bool("geog_to_cart", default=True)
        geogToCart.meta['tip'] = "Convert geographic to Cartesian."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="geo2cart"):
        Application.__init__(self, name)
        self.numPoles = 0
        self.polesSource = []
        return

    def main(self):
        # import pdb
        # pdb.set_trace()
        self._readPoles()
        f = open(self.polesOutputFile, 'w')
        if self.geogToCart:
            self._writeCartPoles(f)
        else:
            self._writeGeogPoles(f)
        f.close()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        self.polesInputFile = self.inventory.polesInputFile
        self.polesOutputFile = self.inventory.polesOutputFile
        self.geogToCart = self.inventory.geogToCart
        self.spaceDim = 3
        return

    def _writeCartPoles(self, f):
        """Computes Cartesian poles and writes to output file.
        """

        iCount = 0
        for pole in range(self.numPoles):
            lat = math.radians(self.polesSource[iCount])
            lon = math.radians(self.polesSource[iCount+1])
            omega = self.polesSource[iCount+2]
            omegax = omega * math.cos(lat) * math.cos(lon)
            omegay = omega * math.cos(lat) * math.sin(lon)
            omegaz = omega * math.sin(lat)
            f.write(' %.12e' % omegax)
            f.write(' %.12e' % omegay)
            f.write(' %.12e' % omegaz)
            f.write('\n')
            iCount += 3
        return

    def _writeGeogPoles(self, f):
        """Computes geographic poles and writes to output file.
        """

        iCount = 0
        for pole in range(self.numPoles):
            omegax = self.polesSource[iCount]
            omegay = self.polesSource[iCount+1]
            omegaz = self.polesSource[iCount+2]
            lat = math.degrees(math.atan2(
                omegaz, math.sqrt(omegax*omegax + omegay*omegay)))
            lon = math.degrees(math.atan2(omegay, omegax))
            omega = math.sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz)
            f.write(' %.12e' % lat)
            f.write(' %.12e' % lon)
            f.write(' %.12e' % omega)
            f.write('\n')
            iCount += 3
        return

    def _readPoles(self):
        """Reads poles from a file.
        """
        f = file(self.polesInputFile)
        for line in f.readlines():
            if not line.startswith('#'):
                data = line.split()
                for dim in range(self.spaceDim):
                    self.polesSource.append(float(data[dim]))
                self.numPoles += 1
        f.close()
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = PoleConvert()
    app.run()

# End of file
