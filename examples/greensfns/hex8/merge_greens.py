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

## @file greensfns/merge_greens

## @brief Python application to merge Green's functions that were created using
## separate fault patches with gfgen.py (and then PyLith). Functions with
## duplicate coordinates are removed and all the necessary files are renumbered
## and moved to a specified location. A metadata file is produced, along with
## a file showing the correspondence between the original and moved files.

import math
import numpy
import os
import re
import glob
from pyre.units.time import s
from pyre.units.length import m

from pyre.applications.Script import Script as Application

class MergeGreens(Application):
  """
  Python application to merge Green's functions that were created using
  separate fault patches with gfgen.py (and then PyLith). Functions with
  duplicate coordinates are removed and all the necessary files are renumbered
  and moved to a specified location. A metadata file is produced, along with
  a file showing the correspondence between the original and moved files.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing MergeGreens facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MergeGreens facilities and properties.
    ##
    ## \b Properties
    ## @li \b metadata_input File with a list of all the metadata input files.
    ## @li \b impulse_input File with root names of all input impulse files.
    ## @li \b response_input File with root names of all input response files.
    ## @li \b metadata_output Name of output metadata file.
    ## @li \b correspondence_output Name of output correspondence file.
    ## @li \b impulse_output Root name of impulse output files.
    ## @li \b response_output Root name of response output files.
    ## @li \b impulse_input_width Width of input impulse number field.
    ## @li \b impulse_output_width Width of output impulse number field.
    ## @li \b distance_tol Distance tolerance to determine coincident vertices.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    metadataInput = pyre.inventory.str("metadata_input",
                                       default="metadata.list")
    metadataInput.meta['tip'] = "File with a list of input metadata files."

    impulseInput = pyre.inventory.str("impulse_input",
                                      default="impulse.list")
    impulseInput.meta['tip'] = "File with a list of input impulse root names."

    responseInput = pyre.inventory.str("response_input",
                                       default="response.list")
    responseInput.meta['tip'] = "File with a list of input response root names."

    metadataOutput = pyre.inventory.str("metadata_output",
                                        default="impulse_description.txt")
    metadataOutput.meta['tip'] = "Name of output file with merged metadata."

    correspondenceOutput = pyre.inventory.str("correspondence_output",
                                              default="correspondence.txt")
    correspondenceOutput.meta['tip'] = "Name of output correspondence file."

    impulseOutput = pyre.inventory.str("impulse_output",
                                       default="gfimpulse.vtk")
    impulseOutput.meta['tip'] = "Root name of impulse output files."

    responseOutput = pyre.inventory.str("response_output",
                                        default="gfresponse.vtk")
    responseOutput.meta['tip'] = "Root name of response output files."
    
    impulseInputWidth = pyre.inventory.int("impulse_input_width", default=4)
    impulseInputWidth.meta['tip'] = "Width of input impulse number field."
    
    impulseOutputWidth = pyre.inventory.int("impulse_output_width", default=5)
    impulseOutputWidth.meta['tip'] = "Width of output impulse number field."
    
    distanceTol = pyre.inventory.dimensional("distance_tol", default=0.1*m)
    distanceTol.meta['tip'] = "Distance tolerance for coincident vertices."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="merge_greens"):
    Application.__init__(self, name)

    self.numPatches = 0
    self.numOriginalImpulses = 0
    self.numMergedImpulses = 0
    self.numOutputImpulses = 0

    self.impulseOutputDir = ''
    self.responseOutputDir = ''
    self.impulseOutputRoot = ''
    self.responseOutputRoot = ''

    self.numImpulsesPerPatch = []
    self.metadataInputFiles = []
    self.impulseInputRoots = []
    self.responseInputRoots = []
    self.metadata = []

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readLists()
    self._setupOutput()
    self._readMetadata()
    self._excludeImpulses()

    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)

    # File info.
    self.metadataInput = self.inventory.metadataInput
    self.impulseInput = self.inventory.impulseInput
    self.responseInput = self.inventory.responseInput
    self.metadataOutput = self.inventory.metadataOutput
    self.correspondenceOutput = self.inventory.correspondenceOutput
    self.impulseOutput = self.inventory.impulseOutput
    self.responseOutput = self.inventory.responseOutput
    
    # Impulse information
    self.impulseInputWidth = self.inventory.impulseInputWidth
    self.impulseOutputWidth = self.inventory.impulseOutputWidth
    
    # Excluded vertex information
    self.distanceTol = self.inventory.distanceTol.value
    
    return
      

  def _readLists(self):
    """
    Function to read lists for input files.
    """
    import os
    import re
    
    # Read names of metadata files.
    f = open(self.metadataInput, 'r')
    lines = f.readlines()
    for line in lines:
      self.metadataInputFiles.append(line.rstrip('\n'))
      self.numPatches += 1
    f.close()

    # Read root filenames of impulse files.
    f = open(self.impulseInput, 'r')
    lines = f.readlines()
    for line in lines:
      impulseRoot = line.rstrip('\n')
      totalInputPath = os.path.normpath(os.path.join(os.getcwd(), impulseRoot))
      inputDir = os.path.dirname(totalInputPath)
      baseInputName = os.path.basename(totalInputPath)
      baseInputNameLen = len(baseInputName)
      baseInputNameSuffStripped = baseInputName
      if baseInputName.endswith(".vtk"):
        baseInputNameSuffStripped = baseInputName[0:baseInputNameLen-4]
      baseInputNameTimeStripped = baseInputNameSuffStripped
      testFind = re.search('_t[0-9]*$', baseInputNameSuffStripped)
      if testFind != None:
        timeInd = baseInputNameSuffStripped.rfind(testFind.group(0))
        baseInputNameTimeStripped = baseInputNameSuffStripped[0:timeInd]
      self.impulseInputRoots.append(os.path.join(inputDir,
                                                 baseInputNameTimeStripped))
    f.close()

    # Read root filenames of response files.
    f = open(self.responseInput, 'r')
    lines = f.readlines()
    for line in lines:
      responseRoot = line.rstrip('\n')
      totalInputPath = os.path.normpath(os.path.join(os.getcwd(), responseRoot))
      inputDir = os.path.dirname(totalInputPath)
      baseInputName = os.path.basename(totalInputPath)
      baseInputNameLen = len(baseInputName)
      baseInputNameSuffStripped = baseInputName
      if baseInputName.endswith(".vtk"):
        baseInputNameSuffStripped = baseInputName[0:baseInputNameLen-4]
      baseInputNameTimeStripped = baseInputNameSuffStripped
      testFind = re.search('_t[0-9]*$', baseInputNameSuffStripped)
      if testFind != None:
        timeInd = baseInputNameSuffStripped.rfind(testFind.group(0))
        baseInputNameTimeStripped = baseInputNameSuffStripped[0:timeInd]
      self.responseInputRoots.append(os.path.join(inputDir,
                                                 baseInputNameTimeStripped))
    f.close()

    print "Number of impulse patches:  %i" % self.numPatches

    return


  def _setupOutput(self):
    """
    Function to setup information for moving files.
    """
    import os
    import re
    
    # Setup root filenames for output impulse files.
    totalOutputPath = os.path.normpath(os.path.join(os.getcwd(),
                                                    self.impulseOutput))
    self.impulseOutputDir = os.path.dirname(totalOutputPath)
    baseOutputName = os.path.basename(totalOutputPath)
    baseOutputNameLen = len(baseOutputName)
    baseOutputNameSuffStripped = baseOutputName
    if baseOutputName.endswith(".vtk"):
      baseOutputNameSuffStripped = baseOutputName[0:baseOutputNameLen-4]
    baseOutputNameTimeStripped = baseOutputNameSuffStripped
    testFind = re.search('_t[0-9]*$', baseOutputNameSuffStripped)
    if testFind != None:
      timeInd = baseOutputNameSuffStripped.rfind(testFind.group(0))
      baseOutputNameTimeStripped = baseOutputNameSuffStripped[0:timeInd]
    self.impulseOutputRoot = os.path.join(self.impulseOutputDir,
                                          baseOutputNameTimeStripped)
    
    # Setup root filenames for output response files.
    totalOutputPath = os.path.normpath(os.path.join(os.getcwd(),
                                                    self.responseOutput))
    self.responseOutputDir = os.path.dirname(totalOutputPath)
    baseOutputName = os.path.basename(totalOutputPath)
    baseOutputNameLen = len(baseOutputName)
    baseOutputNameSuffStripped = baseOutputName
    if baseOutputName.endswith(".vtk"):
      baseOutputNameSuffStripped = baseOutputName[0:baseOutputNameLen-4]
    baseOutputNameTimeStripped = baseOutputNameSuffStripped
    testFind = re.search('_t[0-9]*$', baseOutputNameSuffStripped)
    if testFind != None:
      timeInd = baseOutputNameSuffStripped.rfind(testFind.group(0))
      baseOutputNameTimeStripped = baseOutputNameSuffStripped[0:timeInd]
    self.responseOutputRoot = os.path.join(self.responseOutputDir,
                                           baseOutputNameTimeStripped)

    return
      

  def _readMetadata(self):
    """
    Function to read metadata from each metadata file.
    Results are stored as a list of numpy arrays.
    """
    fileNum = 0
    for metadataFile in self.metadataInputFiles:
      self.metadata.append(numpy.loadtxt(metadataFile, dtype=numpy.float64,
                                         skiprows=1))
      self.numImpulsesPerPatch.append(self.metadata[fileNum].shape[0])
      self.numOriginalImpulses += self.numImpulsesPerPatch[fileNum]
      print "Number of impulses in patch %i:  %i" % (fileNum,
                                                     self.numImpulsesPerPatch[fileNum])
      fileNum += 1

    print "Total number of original impulses:  %i" % self.numOriginalImpulses

    return


  def _getFilename(self, fileRoot, impulse, impulseWidth):
    """
    Function to create a filename given the root filename,
    the impulse number, and the timestamp width.
    """
    impulseNum = int(impulse)
    impulseString = repr(impulseNum).rjust(impulseWidth, '0')
    filename = fileRoot + "_t" + impulseString + ".vtk"
    
    return filename

  def _excludeImpulses(self):
    """
    Function to remove redundant impulses, create the metadata and
    correspondence files, and move files to the proper locations.
    """
    import scipy.spatial.distance
    import os
    
    newline = '\n'
    tab = '\t'
    impulseInFmt = '%0' + str(self.impulseInputWidth) + 'i'
    impulseOutFmt = '%0' + str(self.impulseOutputWidth) + 'i'
    cFmt = '   ' + impulseInFmt + '   ' + impulseOutFmt + newline
    mFmt = impulseOutFmt + 12 * '\t%e' + newline

    # Create output directories for impulses and responses if they don't exist.
    if not os.path.isdir(self.impulseOutputDir):
      os.makedirs(self.impulseOutputDir)
    if not os.path.isdir(self.responseOutputDir):
      os.makedirs(self.responseOutputDir)

    # The first patch retains all impulses.
    metadata = self.metadata[0]
    c = open(self.correspondenceOutput, 'w')
    m = open(self.metadataOutput, 'w')
    patch = 0
    self.numOutputImpulses = self.numImpulsesPerPatch[patch]
    cPatch = 'Patch ' + str(patch) + newline
    cHead = 'Old ID    New ID' + newline
    mHead = "Impulse #" + tab + "X-Coord" + tab + "Y-Coord" + tab + \
            "Z-Coord" + tab + "Normal-X" + tab + "Normal-Y" + tab + \
            "Normal-Z" + tab + "Strike-X" + tab + "Strike-Y" + tab + \
            "Strike-Z" + tab + "Dip-X" + tab + "Dip-Y" + tab + "Dip-Z" \
            + newline
    c.write(cPatch)
    c.write(cHead)
    m.write(mHead)
    newId = 0
    patchId = 0
    for impulse in range(self.numImpulsesPerPatch[patch]):
      cOut = cFmt % (patchId, newId)
      c.write(cOut)
      mOut = mFmt % (int(round(metadata[impulse,0])),
                     metadata[impulse, 1], metadata[impulse, 2],
                     metadata[impulse, 3], metadata[impulse, 4],
                     metadata[impulse, 5], metadata[impulse, 6],
                     metadata[impulse, 7], metadata[impulse, 8],
                     metadata[impulse, 9], metadata[impulse,10],
                     metadata[impulse,11], metadata[impulse,12])
      m.write(mOut)
      impulseFrom = self._getFilename(self.impulseInputRoots[patch],
                                      patchId, self.impulseInputWidth)
      impulseTo = self._getFilename(self.impulseOutputRoot,
                                    newId, self.impulseOutputWidth)
      os.rename(impulseFrom, impulseTo)
      responseFrom = self._getFilename(self.responseInputRoots[patch],
                                       patchId, self.impulseInputWidth)
      responseTo = self._getFilename(self.responseOutputRoot,
                                     newId, self.impulseOutputWidth)
      os.rename(responseFrom, responseTo)
      newId += 1
      patchId += 1

    # Loop over remaining patches.
    for patch in range(1, self.numPatches):
      patchId = 0
      coordsCurrent = metadata[:,1:4]
      coordsPatch = self.metadata[patch][:,1:4]
      distance = scipy.spatial.distance.cdist(coordsPatch, coordsCurrent,
                                              'euclidean')
      minIndices = numpy.argmin(distance, axis=1)
      cPatch = 'Patch ' + str(patch) + newline
      c.write(cPatch)
      c.write(cHead)
      for impulse in range(self.numImpulsesPerPatch[patch]):
        newOutputId = -minIndices[impulse]
        if(distance[impulse, minIndices[impulse]] > self.distanceTol):
          newOutputId = newId
          metadata = numpy.vstack((metadata, self.metadata[patch][impulse,:]))
          mOut = mFmt % (newId,
                         metadata[newId, 1], metadata[newId, 2],
                         metadata[newId, 3], metadata[newId, 4],
                         metadata[newId, 5], metadata[newId, 6],
                         metadata[newId, 7], metadata[newId, 8],
                         metadata[newId, 9], metadata[newId,10],
                         metadata[newId,11], metadata[newId,12])
          m.write(mOut)
          impulseFrom = self._getFilename(self.impulseInputRoots[patch],
                                          patchId, self.impulseInputWidth)
          impulseTo = self._getFilename(self.impulseOutputRoot,
                                        newId, self.impulseOutputWidth)
          os.rename(impulseFrom, impulseTo)
          responseFrom = self._getFilename(self.responseInputRoots[patch],
                                           patchId, self.impulseInputWidth)
          responseTo = self._getFilename(self.responseOutputRoot,
                                         newId, self.impulseOutputWidth)
          os.rename(responseFrom, responseTo)
          newId += 1
          self.numOutputImpulses += 1
        else:
          self.numMergedImpulses += 1
        cOut = cFmt % (patchId, newOutputId)
        c.write(cOut)
        patchId += 1

    print 'Number of output impulses: %i' % self.numOutputImpulses
    print 'Number of merged impulses: %i' % self.numMergedImpulses
    c.close()
    m.close()

    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = MergeGreens()
  app.run()

# End of file
