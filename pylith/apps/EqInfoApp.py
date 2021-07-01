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
# @file pylith/apps/EqInfoApp.py
#
# @brief Python PyLith application

from pythia.pyre.applications.Script import Script as Application

import math
import os

import h5py
import numpy


# ======================================================================
class RuptureStats(object):
    """Python object to hold rupture stats.
    """

    def __init__(self, nsnapshots):
        self.fault = None
        self.timestamp = numpy.zeros((nsnapshots,), dtype=numpy.float64)
        self.ruparea = numpy.zeros((nsnapshots,), dtype=numpy.float64)
        self.potency = numpy.zeros((nsnapshots,), dtype=numpy.float64)
        self.moment = numpy.zeros((nsnapshots,), dtype=numpy.float64)
        return

    def update(self, isnapshot, timestamp, ruparea, potency, moment):
        self.timestamp[isnapshot] = timestamp
        self.ruparea[isnapshot] = ruparea
        self.potency[isnapshot] = potency
        self.moment[isnapshot] = moment
        self.recalculate()
        return

    def recalculate(self):
        self.avgslip = self.potency / (self.ruparea + 1.0e-30)
        self.mommag = -1.0e+30 * numpy.ones(self.moment.shape)
        mask = self.moment > 0.0
        self.mommag[mask] = 2.0 / 3.0 * (numpy.log10(self.moment[mask]) - 9.05)
        return

    def writeObj(self, fout):
        fout.write("class RuptureStats(object):\n"
                   "    pass\n")
        return

    def write(self, fout):
        fout.write("%s = RuptureStats()\n" % self.fault)

        self._writeArray("timestamp", fout)
        self._writeArray("ruparea", fout)
        self._writeArray("potency", fout)
        self._writeArray("moment", fout)
        self._writeArray("avgslip", fout)
        self._writeArray("mommag", fout)

        return

    def _writeArray(self, name, fout):
        vals = self.__getattribute__(name)
        for i in range(len(vals)):
            if math.isnan(vals[i]) or math.isinf(vals[i]):
                if vals[i] > 0:
                    vals[i] = 1.0e+30
                else:
                    vals[i] = -1.0e+30
        g = ("%14.6e" % v for v in vals)
        astr = ", ".join(g)
        fout.write("%s.%s = [%s]\n" % (self.fault, name, astr))
        return


# ======================================================================
# EqInfoApp class
class EqInfoApp(Application):
    """Python EqInfoApp application.
    """

    import pythia.pyre.inventory

    faults = pythia.pyre.inventory.list("faults", default=[])
    faults.meta['tip'] = "Array of fault names."

    filenamePattern = pythia.pyre.inventory.str(
        "filename_pattern", default="output/fault_%s.h5")
    filenamePattern.meta['tip'] = "Pattern for fault files."

    snapshots = pythia.pyre.inventory.list("snapshots", default=[-1])
    snapshots.meta['tip'] = "Array of timestamps for slip snapshots (-1 == last time step)."

    from pythia.pyre.units.time import second
    snapshotUnits = pythia.pyre.inventory.dimensional(
        "snapshot_units", default=1 * second)
    snapshotUnits.meta['tip'] = "Units for timestamps in array of snapshots."

    filenameOut = pythia.pyre.inventory.str(
        "output_filename", default="eqstats.py")
    filenameOut.meta['tip'] = "Filename for output."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbProps = pythia.pyre.inventory.facility(
        "db_properties", family="spatial_database", factory=SimpleDB)
    dbProps.meta['tip'] = "Spatial database for elastic properties."

    from spatialdata.geocoords.CSCart import CSCart
    cs = pythia.pyre.inventory.facility(
        "coordsys", family="coordsys", factory=CSCart)
    cs.meta['tip'] = "Coordinate system associated with mesh."

    typos = pythia.pyre.inventory.str("typos", default="pedantic",
                                      validator=pythia.pyre.inventory.choice(['relaxed', 'strict', 'pedantic']))
    typos.meta['tip'] = "Specifies the handling of unknown properties and " \
        "facilities"

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="eqinfoapp"):
        """Constructor.
        """
        Application.__init__(self, name)
        return

    def main(self, *args, **kwds):
        """Run the application.
        """
        nfaults = len(self.faults)

        if nfaults == 0:
            raise ValueError("No faults specified")

        nsnapshots = len(self.snapshots)
        statsFaults = [None] * nfaults

        ifault = 0
        for fault in self.faults:
            filenameIn = self.filenamePattern % fault
            if not os.path.isfile(filenameIn):
                raise IOError(
                    "Could not open PyLith fault output file '%s'." % filenameIn)

            h5 = h5py.File(filenameIn, "r", driver='sec2')
            vertices = h5['geometry/vertices'][:]
            cells = numpy.array(h5['topology/cells'][:], dtype=numpy.int32)
            timestamps = h5['time'][:]
            cellsArea = self._calcCellArea(cells, vertices)
            cellsShearMod = self._getShearModulus(cells, vertices)

            stats = RuptureStats(nsnapshots)
            statsFaults[ifault] = stats
            stats.fault = fault

            isnapshot = 0
            for snapshot in self.snapshots:
                # Get slip at snapshot
                istep = self._findTimeStep(snapshot, timestamps)
                slip = h5['vertex_fields/slip'][istep, :, :]
                if len(slip.shape) > 2:
                    slip = slip.squeeze(axis=0)

                cellsSlip = self._ptsToCells(slip, cells)
                cellsSlipMag = self._vectorMag(cellsSlip)
                mask = cellsSlipMag > 0.0

                ruparea = numpy.sum(cellsArea[mask])
                potency = numpy.sum(cellsSlipMag * cellsArea)
                moment = numpy.sum(cellsSlipMag * cellsArea * cellsShearMod)

                stats.update(
                    isnapshot, timestamp=timestamps[istep], ruparea=ruparea, potency=potency, moment=moment)

                isnapshot += 1
            h5.close()
            ifault += 1

        statsTotal = RuptureStats(nsnapshots)
        statsTotal.fault = "all"
        isnapshot = 0
        for snapshot in self.snapshots:
            istep = self._findTimeStep(snapshot, timestamps)
            statsTotal.timestamp[isnapshot] = timestamps[istep]
            isnapshot += 1

        ruparea = statsTotal.ruparea
        potency = statsTotal.potency
        moment = statsTotal.moment
        for s in statsFaults:
            ruparea += s.ruparea
            potency += s.potency
            moment += s.moment

        statsTotal.recalculate()

        fout = open(self.filenameOut, "w")
        statsTotal.writeObj(fout)
        statsTotal.write(fout)
        for s in statsFaults:
            s.write(fout)
        fout.close()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        self.snapshots = list(map(float, self.snapshots))

        return

    def _calcCellArea(self, cells, vertices):
        (ncells, ncorners) = cells.shape
        if ncorners == 1:  # point
            area = numpy.ones((ncells,), dtype=numpy.float64)
        elif ncorners == 2:  # line2
            area = self._vectorMag(
                vertices[cells[:, 1]] - vertices[cells[:, 0]])
        elif ncorners == 3:  # tri3
            v01 = vertices[cells[:, 1]] - vertices[cells[:, 0]]
            v02 = vertices[cells[:, 2]] - vertices[cells[:, 0]]
            area = 0.5 * self._vectorMag(numpy.cross(v01, v02))
        elif ncorners == 4:  # quad4
            v01 = vertices[cells[:, 1]] - vertices[cells[:, 0]]
            v02 = vertices[cells[:, 2]] - vertices[cells[:, 0]]
            v03 = vertices[cells[:, 3]] - vertices[cells[:, 0]]
            area = 0.5 * self._vectorMag(numpy.cross(v01, v02)) + \
                0.5 * self._vectorMag(numpy.cross(v02, v03))
        else:
            raise ValueError(
                "Unknown case for number of cell corners (%d)." % ncorners)
        return area

    def _getShearModulus(self, cells, vertices):
        coords = self._ptsToCells(vertices, cells)
        db = self.dbProps
        db.open()
        db.setQueryValues(["density", "vs"])
        (ncells, ndims) = coords.shape
        data = numpy.zeros((ncells, 2), dtype=numpy.float64)
        err = numpy.zeros((ncells,), dtype=numpy.int32)
        db.multiquery(data, err, coords, self.cs)
        db.close()
        shearMod = data[:, 0] * data[:, 1]**2
        return shearMod

    def _findTimeStep(self, value, timestamps):
        if value == -1:
            i = len(timestamps) - 1
        else:
            tdiff = numpy.abs(timestamps - value * self.snapshotUnits.value)
            mindiff = numpy.min(tdiff)
            i = numpy.where(tdiff < mindiff + 1.0e-10)[0]
            if len(i) > 1:
                raise ValueError(
                    "Could not find snapshot %12.4e s in time stamps." % value)
        return i

    def _vectorMag(self, v):
        (npts, ndims) = v.shape
        mag = numpy.zeros((npts,), dtype=numpy.float64)
        for i in range(ndims):
            mag += v[:, i]**2
        mag = mag**0.5
        return mag

    def _ptsToCells(self, valueP, cells):
        (ncells, ncorners) = cells.shape
        if len(valueP.shape) > 1:
            (nvertices, nvals) = valueP.shape
            valueC = numpy.zeros((ncells, nvals), dtype=numpy.float64)
            for i in range(ncorners):
                valueC[:, :] += valueP[cells[:, i], :]
        else:
            nvertices = valueP.shape
            valueC = numpy.zeros((ncells,), dtype=numpy.float64)
            for i in range(ncorners):
                valueC[:] += valueP[cells[:, i]]
        valueC /= ncorners
        return valueC


# End of file
