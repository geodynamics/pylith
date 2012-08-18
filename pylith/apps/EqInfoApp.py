#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/apps/PyLithApp.py
##
## @brief Python PyLith application

from pyre.applications.Script import Script as Application

import h5py
import numpy
import os

# ======================================================================
class RuptureStats(object):
  """
  Python object to hold rupture stats.
  """

  def __init__(self, nsnapshots):
    self.fault = None
    self.timestamp = numpy.zeros( (nsnapshots,), dtype=numpy.float64)
    self.ruparea = numpy.zeros( (nsnapshots,), dtype=numpy.float64)
    self.potency = numpy.zeros( (nsnapshots,), dtype=numpy.float64)
    self.moment = numpy.zeros( (nsnapshots,), dtype=numpy.float64)
    return

  def update(self, isnapshot, timestamp, ruparea, potency, moment):
    self.timestamp[isnapshot] = timestamp
    self.ruparea[isnapshot] = ruparea
    self.potency[isnapshot] = potency
    self.moment[isnapshot] = moment
    self.recalculate()
    return


  def recalculate(self):
    self.avgslip = self.potency / self.ruparea
    self.mommag = 2.0/3.0*numpy.log10(self.moment) - 10.7
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

    g = ("%14.6e" % v for v in self.__getattribute__(name))
    astr = ", ".join(g)
    fout.write("%s.%s = [%s]\n" % (self.fault, name, astr))
    return


# ======================================================================
# EqInfoApp class
class EqInfoApp(Application):
  """
  Python EqInfoApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Application.Inventory):
    """
    Python object for managing EqInfoApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing EqInfoApp facilities and properties.
    ##
    ## \b Properties
    ## @li \b faults Array of fault names.
    ## @li \b filename_pattern Pattern for fault files.
    ## @li \b snapshots Array of timestamps for slip snapshots.
    ## @li \b snapshotUnits Units for timestamps in array of snapshots.
    ## @li \b output_filename Filename for output.
    ##
    ## \b Facilities
    ## @li \b db_properties Spatial database for elastic properties.
    ## @li \b coordsys Coordinate system associated with mesh.

    import pyre.inventory

    faults = pyre.inventory.list("faults", default=[])
    faults.meta['tip'] = "Array of fault names."

    filenamePattern = pyre.inventory.str("filename_pattern", default="output/fault_%s.h5")
    filenamePattern.meta['tip'] = "Pattern for fault files."

    snapshots = pyre.inventory.list("snapshots", default=[-1])
    snapshots.meta['tip'] = "Array of timestamps for slip snapshots (-1 == last time step)."

    from pyre.units.time import second
    snapshotUnits = pyre.inventory.dimensional("snapshot_units", default=1*second)
    snapshotUnits.meta['tip'] = "Units for timestamps in array of snapshots."
    
    filenameOut = pyre.inventory.str("output_filename", default="eqstats.py")
    filenameOut.meta['tip'] = "Filename for output."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbProps = pyre.inventory.facility("db_properties", family="spatial_database", factory=SimpleDB)
    dbProps.meta['tip'] = "Spatial database for elastic properties."
    
    from spatialdata.geocoords.CSCart import CSCart
    cs = pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    cs.meta['tip'] = "Coordinate system associated with mesh."
    
    typos = pyre.inventory.str("typos", default="pedantic",
                               validator=pyre.inventory.choice(['relaxed', 'strict', 'pedantic']))
    typos.meta['tip'] = "Specifies the handling of unknown properties and " \
        "facilities"
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="eqinfoapp"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    return


  def main(self, *args, **kwds):
    """
    Run the application.
    """
    nfaults = len(self.faults)
    
    if nfaults == 0:
      raise ValueError("No faults specified")

    nsnapshots = len(self.snapshots)
    statsFaults = [RuptureStats(nsnapshots)]*nfaults
    
    ifault = 0
    for fault in self.faults:
      filenameIn = self.filenamePattern % fault
      if not os.path.isfile(filenameIn):
        raise IOError("Could not open PyLith fault output file '%s'." % filenameIn)

      h5 = h5py.File(filenameIn, "r", driver='sec2')
      vertices = h5['geometry/vertices'][:]
      cells = h5['topology/cells'][:]
      timestamps = h5['time'][:]
      cellsArea = self._calcCellArea(cells, vertices)
      cellsShearMod = self._getShearModulus(cells, vertices)

      stats = statsFaults[ifault]
      stats.fault = fault
      
      isnapshot = 0
      for snapshot in self.snapshots:
        # Get slip at snapshot
        istep = self._findTimeStep(snapshot, timestamps)
        slip = h5['vertex_fields/slip'][istep,:,:]
        
        slipMag = self._vectorMag(slip)
        cellsSlipMag = self._ptsToCells(slipMag, cells)
        mask = cellsSlipMag > 0.0

        ruparea = numpy.sum(cellsArea[mask])
        potency = numpy.sum(cellsSlipMag*cellsArea)
        moment = numpy.sum(cellsSlipMag*cellsArea*cellsShearMod)

        stats.update(isnapshot, timestamp=snapshot, ruparea=ruparea, potency=potency, moment=moment)
        
        isnapshot += 1
      h5.close()
      ifault += 1
      
    statsTotal = RuptureStats(nsnapshots)
    statsTotal.fault = "all"
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
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.faults = self.inventory.faults
    self.filenamePattern = self.inventory.filenamePattern
    self.snapshots = map(float, self.inventory.snapshots)
    self.snapshotUnits = self.inventory.snapshotUnits
    self.dbProps = self.inventory.dbProps
    self.cs = self.inventory.cs
    self.filenameOut = self.inventory.filenameOut
    self.typos = self.inventory.typos

    return


  def _calcCellArea(self, cells, vertices):
    (ncells, ncorners) = cells.shape
    if ncorners == 1:
      area = numpy.ones( (ncells,), dtype=numpy.float64)
    elif ncorners == 2:
      area = self._vectorMag(vertices[cells[:,1]] - vertices[cells[:,0]])
    elif ncorners == 3:
      v01 = vertices[cells[:,1]] - vertices[cells[:,0]]
      v02 = vertices[cells[:,2]] - vertices[cells[:,0]]
      area = 0.5*self._vectorMag(numpy.cross(v01, v02))
    elif ncorners == 4:
      v01 = vertices[cells[:,1]] - vertices[cells[:,0]]
      v02 = vertices[cells[:,2]] - vertices[cells[:,0]]
      area = self._vectorMag(numpy.cross(v01, v02))
    else:
      raise ValueError("Unknown case for number of cell corners (%d)." % ncorners)
    return area


  def _getShearModulus(self, cells, vertices):
    coords = self._ptsToCells(vertices, cells)
    db = self.dbProps
    db.open()
    db.queryVals(["density","vs"])
    (ncells, ndims) = coords.shape
    data = numpy.zeros( (ncells, 2), dtype=numpy.float64)
    err = numpy.zeros( (ncells,), dtype=numpy.int32)
    db.multiquery(data, err, coords, self.cs)
    db.close()
    shearMod = data[:,0]*data[:,1]**2
    return shearMod


  def _findTimeStep(self, value, timestamps):
    if value == -1:
      i = len(timestamps)-1
    else:
      tdiff = numpy.abs(timestamps-value*self.snapshotUnits.value)
      mindiff = numpy.min(tdiff)
      i = numpy.where(tdiff < mindiff+1.0e-10)[0]
      if len(i) > 1:
        raise ValueError("Could not find snapshot %12.4e s in time stamps." % value)
    return i


  def _vectorMag(self, v):
    (npts, ndims) = v.shape
    mag = numpy.zeros( (npts,), dtype=numpy.float64)
    for i in xrange(ndims):
      mag += v[:,i]**2
    mag = mag**0.5
    return mag


  def _ptsToCells(self, valueP, cells):
    (ncells, ncorners) = cells.shape
    if len(valueP.shape) > 1:
      (nvertices, nvals) = valueP.shape
      valueC = numpy.zeros( (ncells,nvals), dtype=numpy.float64)
      for i in xrange(ncorners):
        valueC[:,:] += valueP[cells[:,i],:]
    else:
      nvertices = valueP.shape
      valueC = numpy.zeros( (ncells,), dtype=numpy.float64)
      for i in xrange(ncorners):
        valueC[:] += valueP[cells[:,i]]
    valueC /= ncorners
    return valueC


# End of file 
