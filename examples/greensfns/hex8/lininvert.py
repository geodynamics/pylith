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

## @file greensfns/lininvert

## @brief Python application to perform a linear inversion, given the necessary
## matrices. In addition to performing the inversion, a number of (optional)
## files are produced containing inversion information.

import math
import numpy
import sys
from pyre.units.length import km
from pyre.units.length import m

from pyre.applications.Script import Script as Application

class LinInvert(Application):
  """
  Python application to perform a linear inversion, given the necessary
  matrices. In addition to performing the inversion, a number of (optional)
  files are produced containing inversion information.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing LinInvert facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing LinInvert facilities and properties.
    ##
    ## \b Properties
    ## @li \b data_input_file File containing data values.
    ## @li \b apriori_input_file File containing a priori parameter values.
    ## @li \b data_cov_input_file File containing data covariance.
    ## @li \b apriori_cov_input_file File containing a priori covariance.
    ## @li \b data_design_input_file File containing data design matrix.
    ## @li \b info_output_file Output file with inversion summary.
    ## @li \b data_output_file Output file containing data information.
    ## @li \b param_output_file Output file containing parameter information.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    dataInputFile = pyre.inventory.str("data_input_file", default="data.txt")
    dataInputFile.meta['tip'] = "File containing data values."

    aprioriInputFile = pyre.inventory.str("apriori_input_file",
                                        default="apriori_values.txt")
    aprioriInputFile.meta['tip'] = "File containing a priori parameter values."

    dataCovInputFile = pyre.inventory.str("data_cov_input_file",
                                          default="data_cov.txt")
    dataCovInputFile.meta['tip'] = "File containing data covariance."

    aprioriCovInputFile = pyre.inventory.str("apriori_cov_input_file",
                                             default="apriori_cov.txt")
    aprioriCovInputFile.meta['tip'] = "File containing a priori covariance."

    dataDesignInputFile = pyre.inventory.str("data_design_input_file",
                                             default="data_design.txt")
    dataDesignInputFile.meta['tip'] = "File containing data design matrix."

    infoOutputFile = pyre.inventory.str("info_output_file",
                                        default="inversion_info.txt")
    infoOutputFile.meta['tip'] = "Output file with inversion summary."

    dataOutputFile = pyre.inventory.str("data_output_file",
                                        default="data_output.txt")
    dataOutputFile.meta['tip'] = "Output file containing data information."

    paramOutputFile = pyre.inventory.str("param_output_file",
                                         default="param_output.txt")
    paramOutputFile.meta['tip'] = "Output file with parameter information."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="lininvert"):
    Application.__init__(self, name)

    self.numObs = 0
    self.numParams = 0
    self.chiSquare = 0.0

    self.dataVals = None
    self.dataCovVec = None
    self.dataCov = None
    self.dataDesign = None
    self.aprioriVals = None
    self.aprioriCov = None

    self.hData = None
    self.hApriori = None
    self.resolData = None
    self.solution = None
    self.fData = None
    self.predicted = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readMats()
    self._invertData()
    self._outputDataInfo()
    self._outputParamInfo()
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
    self.dataInputFile = self.inventory.dataInputFile
    self.aprioriInputFile = self.inventory.aprioriInputFile
    self.dataCovInputFile = self.inventory.dataCovInputFile
    self.aprioriCovInputFile = self.inventory.aprioriCovInputFile
    self.dataDesignInputFile = self.inventory.dataDesignInputFile
    self.infoOutputFile = self.inventory.infoOutputFile
    self.dataOutputFile = self.inventory.dataOutputFile
    self.paramOutputFile = self.inventory.paramOutputFile

    return
      

  def _readMats(self):
    """
    Function to read input matrices.
    """
    self.dataVals = numpy.loadtxt(self.dataInputFile, dtype=numpy.float64)
    self.dataCovVec = numpy.loadtxt(self.dataCovInputFile, dtype=numpy.float64)
    self.dataCov = numpy.diag(self.dataCovVec)
    self.dataDesign = numpy.loadtxt(self.dataDesignInputFile,
                                    dtype=numpy.float64)
    self.aprioriVals = numpy.loadtxt(self.aprioriInputFile, dtype=numpy.float64)
    self.aprioriCov = numpy.loadtxt(self.aprioriCovInputFile,
                                    dtype=numpy.float64)

    self.numObs = self.dataVals.shape[0]
    self.numParams = self.aprioriVals.shape[0]
    print "Number of observations: %i" % self.numObs
    print "Number of parameters: %i" % self.numParams
    sys.stdout.flush()

    return


  def _invertData(self):
    """
    Function to perform inversion.
    """
    # Perform inversion.
    covAAdT = numpy.dot(self.aprioriCov, numpy.transpose(self.dataDesign))
    sum1 = self.dataCov + numpy.dot(self.dataDesign, covAAdT)
    self.hData = numpy.dot(covAAdT, numpy.linalg.inv(sum1))
    self.resolData = numpy.dot(self.hData, self.dataDesign)
    self.hApriori = numpy.identity(self.numParams, dtype=numpy.float64) - \
                    self.resolData
    self.solution = numpy.dot(self.hData, self.dataVals) + \
                    numpy.dot(self.hApriori, self.aprioriVals)

    return


  def _outputDataInfo(self):
    """
    Function to write out general info and info related to data.
    """

    # Compute data-related quantities.
    fDataVec = numpy.sqrt(1.0/self.dataCovVec)
    fDataVecInv = 1.0/fDataVec
    self.fData = numpy.diag(fDataVec)
    self.predicted = numpy.dot(self.dataDesign, self.solution)
    dataMisfit = self.dataVals - self.predicted
    dataWeightMisfit = dataMisfit * fDataVec
    self.chiSquare = numpy.sum(dataWeightMisfit * dataWeightMisfit)
    dataMisfitNorm = numpy.linalg.norm(dataMisfit)
    dataWeightMisfitNorm = numpy.linalg.norm(dataWeightMisfit)
    
    # Write out inversion info.
    string1 = "\nNumber of observations:  %i\n" % self.numObs
    string2 = "Number of parameters:  %i\n" % self.numParams
    string3 = "Data Chi-square value:  %e\n" % self.chiSquare
    string4 = "Data residual norm:  %e\n" % dataMisfitNorm
    string5 = "Weighted data residual norm:  %e\n" % dataWeightMisfitNorm
    print string1
    print string2
    print string3
    print string4
    print string5
    sys.stdout.flush()
    i = open(self.infoOutputFile, 'w')
    i.write(string1)
    i.write(string2)
    i.write(string3)
    i.write(string4)
    i.write(string5)
    i.close()
                               
    # Write out data info.
    # Write out the following columns:
    # 1.  Observed data value.
    # 2.  Predicted data value.
    # 3.  Standard deviation of observation.
    # 4.  Observed minus predicted value.
    # 5.  Observed minus predicted weighted by standard deviation.
    
    dataOut = numpy.transpose(numpy.vstack((self.dataVals, self.predicted,
                                            fDataVecInv, dataMisfit,
                                            dataWeightMisfit)))
    numpy.savetxt(self.dataOutputFile, dataOut)

    return


  def _outputParamInfo(self):
    """
    Function to write out parameter-related information.
    """
    
    # Compute a posteriori covariance contributions.
    aposterioriCov = numpy.dot(self.hApriori, self.aprioriCov)
    prod1 = numpy.dot(self.hData, self.dataCov)
    aposterioriCovD = numpy.dot(prod1, numpy.transpose(self.hData))
    aposterioriCovA = numpy.dot(aposterioriCov, numpy.transpose(self.hApriori))

    # Compute eigenvalues and eigenvectors of a priori covariance.
    aprioriEigVals, aprioriEigVecs = numpy.linalg.eigh(self.aprioriCov)
    aprioriValsVec = numpy.sqrt(1.0/aprioriEigVals)
    aprioriValsVecInv = 1.0/aprioriValsVec
    aprioriValsDiag = numpy.diag(aprioriValsVec)
    aprioriValsDiagInv = numpy.diag(aprioriValsVecInv)
    fApriori = numpy.dot(aprioriValsDiag, numpy.transpose(aprioriEigVecs))
    fAprioriInv = numpy.dot(aprioriEigVecs, aprioriValsDiagInv)
    
    # Compute a posteriori covariances in primed coordinate system.
    prod2 = numpy.dot(fApriori, aposterioriCov)
    aposterioriCovPrime = numpy.dot(prod2, numpy.transpose(fApriori))
    hDataPrime = numpy.dot(aposterioriCovPrime,
                           numpy.transpose(self.dataDesign))
    prod3 = numpy.dot(self.dataDesign, fAprioriInv)
    dataDesignPrime = numpy.dot(self.fData, prod3)

    aposterioriCovDPrime = numpy.dot(hDataPrime, numpy.transpose(hDataPrime))
    aposterioriCovAPrime = numpy.dot(aposterioriCovPrime,
                                     numpy.transpose(aposterioriCovPrime))

    resolDataPrime = numpy.dot(hDataPrime, dataDesignPrime)

    # Columns to output:
    # 1.  Predicted parameter value.
    # 2.  A priori parameter value.
    # 3.  A priori covariance.
    # 4.  A posteriori covariance.
    # 5.  A posteriori covariance (normalized).
    # 6.  Contribution of data to a posteriori covariance.
    # 7.  Contribution of data to a posteriori covariance (normalized).
    # 8.  Contribution of a priori data to a posteriori covariance.
    # 9.  Contribution of a priori data to a posteriori covariance (normalized).
    # 10.  Contribution of data to resolution matrix.
    # 11.  Contribution of data to resolution matrix (normalized).

    # Extract necessary diagonals.
    aprioriCovDiag = numpy.diag(self.aprioriCov)
    aposterioriCovDiag = numpy.diag(aposterioriCov)
    aposterioriCovPrimeDiag = numpy.diag(aposterioriCovPrime)
    aposterioriCovDDiag = numpy.diag(aposterioriCovD)
    aposterioriCovDPrimeDiag = numpy.diag(aposterioriCovDPrime)
    aposterioriCovADiag = numpy.diag(aposterioriCovA)
    aposterioriCovAPrimeDiag = numpy.diag(aposterioriCovAPrime)
    resolDataDiag = numpy.diag(self.resolData)
    resolDataPrimeDiag = numpy.diag(resolDataPrime)

    # Output columns.
    paramOut = numpy.transpose(numpy.vstack((self.solution, self.aprioriVals,
                                             aprioriCovDiag, aposterioriCovDiag,
                                             aposterioriCovPrimeDiag,
                                             aposterioriCovDDiag,
                                             aposterioriCovDPrimeDiag,
                                             aposterioriCovADiag,
                                             aposterioriCovAPrimeDiag,
                                             resolDataDiag,
                                             resolDataPrimeDiag)))
                            
    numpy.savetxt(self.paramOutputFile, paramOut)

    return
  

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = LinInvert()
  app.run()

# End of file
