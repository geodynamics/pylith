#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file greensfns/designprior

## @brief Python application to create the a priori arrays necessary for an
## inversion, given the impulse locations, the a priori values, and the
## desired correlation length.
## The matrices produced are the a priori data matrix (vector) and the a priori
## covariance matrix. The a priori design matrix is assumed to be an identity
## matrix with dimensions corresponding to the number of parameters.

import math
import numpy
import sys
from pyre.units.length import km
from pyre.units.length import m

from pyre.applications.Script import Script as Application

class DesignPrior(Application):
  """
  Python application to create the a priori arrays necessary for an
  inversion, given the impulse locations, the a priori values, and the
  desired correlation length.
  The matrices produced are the a priori data matrix (vector) and the a priori
  covariance matrix. The a priori design matrix is assumed to be an identity
  matrix with dimensions corresponding to the number of parameters.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing DesignPrior facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing DesignPrior facilities and properties.
    ##
    ## \b Properties
    ## @li \b metadata_input_file File containing metadata for the inversion.
    ## @li \b data_output_file File containing a priori parameter values.
    ## @li \b cov_output_file Output file for a priori covariance matrix.
    ## @li \b correlation_length Correlation length for covariance.
    ## @li \b std_dev Standard deviation to use for covariance.
    ## @li \b apriori_value A priori parameter values (all set to this value).
    ## @li \b diagonal_frac Additional fraction to add to covariance diagonal.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    metadataInputFile = pyre.inventory.str("metadata_input_file",
                                           default="metadata.txt")
    metadataInputFile.meta['tip'] = "File containing inversion metadata."

    dataOutputFile = pyre.inventory.str("data_output_file",
                                        default="data_vals.txt")
    dataOutputFile.meta['tip'] = "Output file with a priori parameter values."

    covOutputFile = pyre.inventory.str("cov_output_file",
                                       default="data_cov.txt")
    covOutputFile.meta['tip'] = "Output file for a priori covriance matrix."

    correlationLength = pyre.inventory.dimensional("correlation_length",
                                                   default=10.0*km)
    correlationLength.meta['tip'] = "Correlation length for covariance."

    stdDev = pyre.inventory.dimensional("std_dev", default=1.0*m)
    stdDev.meta['tip'] = "Standard deviation for covariance."

    aprioriValue = pyre.inventory.dimensional("apriori_value", default=0.0*m)
    aprioriValue.meta['tip'] = "A priori parameter values."

    diagonalFrac = pyre.inventory.float("diagonal_frac", default=0.1)
    diagonalFrac.meta['tip'] = "Additional fraction to add to covariance diagonal."

    minCovar = pyre.inventory.float("min_covar", default=1.0e-5)
    minCovar.meta['tip'] = "Covariance values less than this are set to zero."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="designprior"):
    Application.__init__(self, name)

    self.numParameters = 0

    self.impulseCoords = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readMetadata()
    self._makeCovar()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    # import pdb
    # pdb.set_trace()

    # File info.
    self.metadataInputFile = self.inventory.metadataInputFile
    self.dataOutputFile = self.inventory.dataOutputFile
    self.covOutputFile = self.inventory.covOutputFile

    # Parameters
    self.correlationLength = self.inventory.correlationLength.value
    self.stdDev = self.inventory.stdDev.value
    self.aprioriValue = self.inventory.aprioriValue.value
    self.diagonalFrac = self.inventory.diagonalFrac
    self.minCovar = self.inventory.minCovar

    return
      

  def _readMetadata(self):
    """
    Function to read information describing parameters and write parameter
    values.
    """
    # Get impulse information.
    f = open(self.metadataInputFile, 'r')
    lines = f.readlines()
    coords = []
    impulseData = False
    for line in lines:
      if impulseData:
        lineSplit = line.split()
        impulseCoord = [float(lineSplit[2]), float(lineSplit[3]),
                        float(lineSplit[4])]
        coords.append(impulseCoord)
        self.numParameters += 1
      testFind = line.find('Column')
      if testFind > -1:
        impulseData = True

    f.close()
    self.impulseCoords = numpy.array(coords, dtype=numpy.float64)

    print "Number of parameters: %i" % self.numParameters

    # Write a priori parameter values.
    parameterVals = self.aprioriValue * numpy.ones(self.numParameters)
    numpy.savetxt(self.dataOutputFile, parameterVals)

    return


  def _distanceFunc(self, distance):
    """
    Function to compute autocorrelation based on distance.
    """
    distanceArg = -0.5 * (distance * distance)/\
                  (self.correlationLength * self.correlationLength)
    covariance = self.stdDev * self.stdDev * math.exp(distanceArg)
    return covariance

  
  def _makeCovar(self):
    """
    Function to create a priori covariance matrix for a given correlation
    length.
    """
    import scipy.spatial.distance
    import scipy.stats.mstats

    # Compute distances between all points.
    # Temporary kludge because the prior approach correlated different types
    # of impulses with each other.
    impulseHalf = self.impulseCoords[0:self.numParameters/2,:]
    distance = scipy.spatial.distance.cdist(impulseHalf,
                                            impulseHalf, 'euclidean')

    # Compute covariance.
    multiplier = -0.5/(self.correlationLength * self.correlationLength)
    distanceArg = multiplier * distance * distance
    # distanceFunc = numpy.vectorize(self._distanceFunc)
    # aprioriCovar = distanceFunc(distance)
    covarMult = self.stdDev * self.stdDev
    aprioriCovarHalf = covarMult * numpy.exp(distanceArg)
    aprioriCovarEmpty = numpy.zeros_like(aprioriCovarHalf)

    # Delete all entries below threshold.
    aprioriCovarHalfThresh = scipy.stats.mstats.threshold(
      aprioriCovarHalf, threshmin=self.minCovar)

    # Join together different pieces.
    aprioriCovarTop = numpy.hstack((aprioriCovarHalfThresh, aprioriCovarEmpty))
    aprioriCovarBot = numpy.hstack((aprioriCovarEmpty, aprioriCovarHalfThresh))
    aprioriCovar = numpy.vstack((aprioriCovarTop, aprioriCovarBot))

    # Add additional amount to diagonal.
    covarDiagAdd = self.diagonalFrac * numpy.diag(aprioriCovar)
    diagIndices = numpy.arange(self.numParameters)
    aprioriCovar[diagIndices, diagIndices] += covarDiagAdd

    # Write covariance to file.
    numpy.savetxt(self.covOutputFile, aprioriCovar)

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = DesignPrior()
  app.run()

# End of file
