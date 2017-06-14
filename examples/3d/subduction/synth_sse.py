#!/usr/bin/env python

## @file synth_sse.py

## @brief Python application to create spatial databases for a synthetic
## SSE with time-varying circular slip.

# import pdb
import math
import numpy
from spatialdata.spatialdb.SimpleGridDB import SimpleGridDB
from spatialdata.spatialdb.SimpleGridAscii import SimpleGridAscii
from spatialdata.spatialdb import TimeHistoryIO
from spatialdata.geocoords.CSGeo import CSGeo
from spatialdata.geocoords.CoordSys import CoordSys
from spatialdata.geocoords.CSGeoProj import CSGeoProj
from spatialdata.geocoords.Converter import convert


from pyre.applications.Script import Script as Application

class SynthSse(Application):
  """
  Python application to create spatial databases for a synthetic
  SSE with time-varying Gaussian slip.
  """
  
  import pyre.inventory
  ## Python object for managing SynthSse facilities and properties.
  ##
  ## \b Properties
  ## @li \b rake Rake of fault slip (degrees).
  ## @li \b slip_center (lon,lat) coordinates of slip center.
  ## @li \b slip_radius Radius of slip region (degrees).
  ## @li \b slip_max Maximum slip value (meters).
  ## @li \b slip_sigma_lon Sigma value for longitude.
  ## @li \b slip_sigma_lat Sigma value for latitude.
  ## @li \b slip_times List of times for which to provide amplitudes.
  ## @li \b slip_time_units Units used for slip times.
  ## @li \b slip_amplitudes List of slip amplitudes.
  ## @li \b grid_lon_range Min and max longitude values for grid.
  ## @li \b grid_lat_range Min and max latitude values for grid.
  ## @li \b grid_incr Grid increment (degrees) for spatial database.
  ## @li \b time_db_filename Name of temporal DB output file.
  ##
  ## \b Facilities
  ## @li \b coordsys Coordinate system for output database.
  ## @li \b iohandler Object for writing spatial database.

  rake = pyre.inventory.float("rake", default=1.0)
  rake.meta['tip'] = "Rake of fault slip (degrees)."

  slipCenter = pyre.inventory.list("slip_center", default=[0.0, 0.0])
  slipCenter.meta['tip'] = "(lon,lat) coordinates of slip center."

  slipRadius = pyre.inventory.float("slip_radius", default=1.0)
  slipRadius.meta['tip'] = "Radius of slip region (degrees)."

  slipMax = pyre.inventory.float("slip_max", default=5.0)
  slipMax.meta['tip'] = "Maximum slip value (meters)."

  slipSigmaLon = pyre.inventory.float("slip_sigma_lon", default=0.2)
  slipSigmaLon.meta['tip'] = "Sigma value for longitude."

  slipSigmaLat = pyre.inventory.float("slip_sigma_lat", default=0.2)
  slipSigmaLat.meta['tip'] = "Sigma value for latitude."

  slipTimes = pyre.inventory.list("slip_times", default=[0.0, 0.5, 1.0])
  slipTimes.meta['tip'] = "List of times for which to provide amplitudes."

  slipTimeUnits = pyre.inventory.str("slip_time_units", default="year")
  slipTimeUnits.meta['tip'] = "Units used for slip times."

  slipAmplitudes = pyre.inventory.list("slip_amplitudes",
                                       default=[0.0, 0.5, 1.0])
  slipAmplitudes.meta['tip'] = "List of slip amplitudes."

  gridLonRange = pyre.inventory.list("grid_lon_range", default=[-123.0, -124.0])
  gridLonRange.meta['tip'] = "Min and max longitude values for grid."

  gridLatRange = pyre.inventory.list("grid_lat_range", default=[45.0, 46.0])
  gridLatRange.meta['tip'] = "Min and max latitude values for grid."

  gridIncr = pyre.inventory.float("grid_incr", default=0.05)
  gridIncr.meta['tip'] = "Sigma value for latitude."

  timeDbFilename = pyre.inventory.str("time_db_filename", default="slip.timedb")
  timeDbFilename.meta['tip'] = "Filename of temporal DB output file."

  coordsys = pyre.inventory.facility("coordsys",
                                     family="coordsys",
                                     factory=CSGeo)
  coordsys.meta['tip'] = "Coordinate system for output database."
  
  iohandler = pyre.inventory.facility("iohandler", family="simplegrid_io",
                                      factory=SimpleGridAscii)
  iohandler.meta['tip'] = "Object for writing spatial database."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="synth_sse"):
    Application.__init__(self, name)

    self.lon = None
    self.lat = None
    self.z = None
    self.grid = None
    self.faultSlip = None

    return


  def main(self):
    # pdb.set_trace()
    self.coordsys._configure()
    self.coordsys.initialize()
                                         
    self._makeGrid()
    self._computeGauss()
    self._writeSpatialdb()
    self._writeTemporaldb()

    return
  

  # PRIVATE METHODS /////////////////////////////////////////////////////


  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)

    return


  def _makeGrid(self):
    """
    Function to create a mesh grid for computations.
    """

    lonMin = float(self.gridLonRange[0])
    lonMax = float(self.gridLonRange[1])
    latMin = float(self.gridLatRange[0])
    latMax = float(self.gridLatRange[1])

    lonDiff = lonMax - lonMin
    latDiff = latMax - latMin
    numLon = int(round(lonDiff/self.gridIncr)) + 1
    numLat = int(round(latDiff/self.gridIncr)) + 1

    self.lon = numpy.linspace(lonMin, lonMax, num=numLon, dtype=numpy.float64)
    self.lat = numpy.linspace(latMin, latMax, num=numLat, dtype=numpy.float64)
    self.z = numpy.zeros(1, dtype=numpy.float64)

    lonGrid, latGrid = numpy.meshgrid(self.lon, self.lat)
    zGrid = numpy.zeros_like(lonGrid)
    self.grid = numpy.column_stack((lonGrid.flatten(), latGrid.flatten(),
                                    zGrid.flatten()))

    return


  def _computeGauss(self):
    """
    Function to compute 2D Gaussian slip distribution.
    """

    lonShift = self.grid[:,0] - float(self.slipCenter[0])
    latShift = self.grid[:,1] - float(self.slipCenter[1])

    distance = numpy.sqrt(lonShift * lonShift + latShift * latShift)
    outside = numpy.where(distance > self.slipRadius)

    lonFac = 0.5 * lonShift * lonShift/(self.slipSigmaLon * self.slipSigmaLon)
    latFac = 0.5 * latShift * latShift/(self.slipSigmaLat * self.slipSigmaLat)

    slip = self.slipMax * numpy.exp(-(lonFac + latFac))
    slip[outside] = 0.0

    rakeRadians = math.radians(self.rake)
    llComp = math.cos(rakeRadians)
    udComp = math.sin(rakeRadians)

    llSlip = llComp * slip
    udSlip = udComp * slip
    opSlip = numpy.zeros_like(llSlip)

    self.faultSlip = numpy.column_stack((llSlip, udSlip, opSlip))
    
    return

  
  def _writeSpatialdb(self):
    """
    Write spatial database with fault slip.
    """

    print "Writing spatial database:"

    llSlipInfo = {'name': "left-lateral-slip",
                  'units': "m",
                  'data': self.faultSlip[:,0]}

    udSlipInfo = {'name': "reverse-slip",
                  'units': "m",
                  'data': self.faultSlip[:,1]}

    openInfo = {'name': "fault-opening",
                  'units': "m",
                  'data': self.faultSlip[:,2]}

    data = {'num-x': self.lon.shape[0],
            'num-y': self.lat.shape[0],
            'num-z': 1,
            'points': self.grid,
            'x': self.lon,
            'y': self.lat,
            'z': self.z,
            'coordsys': self.coordsys,
            'data_dim': 2,
            'values': [llSlipInfo, udSlipInfo, openInfo]}

    self.iohandler.write(data)

    return

  
  def _writeTemporaldb(self):
    """
    Write temporal database with time variation of fault slip.
    """

    time = [float(i) for i in self.slipTimes]
    timeArr = numpy.array(time, dtype=numpy.float64)

    amplitude = [float(i) for i in self.slipAmplitudes]
    amplitudeArr = numpy.array(amplitude, dtype=numpy.float64)

    TimeHistoryIO.write(timeArr, amplitudeArr, self.slipTimeUnits,
                        self.timeDbFilename)

    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = SynthSse()
  app.run()

# End of file
