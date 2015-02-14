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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file savpres_ss/savpres_ss

## @brief Python application to compute the Savage and Prescott [1978]
## solution for an infinitely long strike-slip fault embedded in an
## elastic layer overlying a viscoelastic half-space.

import math
import numpy

from pyre.applications.Script import Script as Application

class Savpres_ss(Application):
  """
  Python application to compute the Savage and Prescott [1978] solution
  for an infinitely long strike-slip fault embedded in an elastic layer
  overlying a viscoelastic half-space.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing Savpres_ss facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Savpres_ss facilities and properties.
    ##
    ## \b Properties
    ## @li \b elas_thick Thickness of elastic layer.
    ## @li \b lock_depth Fault locking depth (<= elas_thick).
    ## @li \b recurrence_time Earthquake recurrence time.
    ## @li \b viscosity Viscosity of viscoelastic half-space.
    ## @li \b shear_modulus Shear modulus of layer and half-space.
    ## @li \b plate_velocity Relative plate velocity (left-lateral) across fault.
    ## @li \b number_cycles Number of earthquake cycles to compute.
    ## @li \b number_steps Number of time steps to compute for each cycle.
    ## @li \b number_terms Number of terms to compute for series solution.
    ## @li \b number_points Number of points at which to compute solution.
    ## @li \b delta_x Distance between computation points.
    ## @li \b x_epsilon Offset for computation point closest to the fault.
    ## @li \b output_displ_vtk Output displacement VTK results?
    ## @li \b output_displ_csv Output displacement CSV results?
    ## @li \b output_vel_vtk Output velocity VTK results?
    ## @li \b output_vel_csv Output velocity CSV results?
    ## @li \b displ_vtk_basename Base name for VTK displacement output files.
    ## @li \b displ_csv_filename Filename for CSV displacement output.
    ## @li \b vel_vtk_basename Base name for VTK velocity output files.
    ## @li \b vel_csv_filename Filename for CSV velocity output.
    ## @li \b displ_scale_factor Scaling factor for output displacements.
    ## @li \b coord_scale_factor Scaling factor for output coordinates.
    ## @li \b vel_scale_factor Scaling factor for output velocities.
    ## @li \b coord_units Units used for output coordinates.
    ## @li \b time_units Time units to use for output filenames.
    ## @li \b time_stamp_width Width of time stamp in output filenames.
    ## @li \b title Title to appear at the top of VTK files.

    import pyre.inventory
    from pyre.units.length import m
    from pyre.units.length import km
    from pyre.units.length import cm
    from pyre.units.time import s
    from pyre.units.time import year
    from pyre.units.pressure import Pa
    from pyre.units.pressure import MPa
    from pyre.units.pressure import GPa

    elasThick = pyre.inventory.dimensional("elas_thick", default=20.0*km)
    elasThick.meta['tip'] = "Thickness of elastic layer."

    lockDepth = pyre.inventory.dimensional("lock_depth", default=10.0*km)
    lockDepth.meta['tip'] = "Fault locking depth (<= elastic thickness)."

    recurrenceTime = pyre.inventory.dimensional("recurrence_time",
                                                default=100.0*year)
    recurrenceTime.meta['tip'] = "Earthquake recurrence time."

    viscosity = pyre.inventory.dimensional("viscosity", default=1.0e18*Pa*s)
    viscosity.meta['tip'] = "Half-space viscosity."

    shearModulus = pyre.inventory.dimensional("shear_modulus", default=30.0*GPa)
    shearModulus.meta['tip'] = "Shear modulus of layer and half-space."

    plateVelocity = pyre.inventory.dimensional("plate_velocity",
                                               default=2.0*cm/year)
    plateVelocity.meta['tip'] = "Relative velocity (left-lateral) across the fault."

    numberCycles = pyre.inventory.int("number_cycles", default=10)
    numberCycles.meta['tip'] = "Number of earthquake cycles."

    numberSteps = pyre.inventory.int("number_steps", default=10)
    numberSteps.meta['tip'] = "Number of steps to compute for each cycle."

    numberTerms = pyre.inventory.int("number_terms", default=20)
    numberTerms.meta['tip'] = "Number of terms to compute for series."

    numberPoints = pyre.inventory.int("number_points", default=100)
    numberPoints.meta['tip'] = "Number of points at which to compute solution."

    deltaX = pyre.inventory.dimensional("delta_x", default=2.0*km)
    deltaX.meta['tip'] = "Distance between computation points."

    xEpsilon = pyre.inventory.dimensional("x_epsilon", default=0.001*m)
    xEpsilon.meta['tip'] = "Offset for computation point closest to the fault."

    outputDisplVTK = pyre.inventory.bool("output_displ_vtk", default=False)
    outputDisplVTK.meta['tip'] = "Output displacement VTK files?"

    outputDisplCSV = pyre.inventory.bool("output_displ_csv", default=True)
    outputDisplCSV.meta['tip'] = "Output displacement CSV file?"

    outputVelVTK = pyre.inventory.bool("output_vel_vtk", default=False)
    outputVelVTK.meta['tip'] = "Output velocity VTK files?"

    outputVelCSV = pyre.inventory.bool("output_vel_csv", default=True)
    outputVelCSV.meta['tip'] = "Output velocity CSV file?"

    displVTKBaseName = pyre.inventory.str("displ_vtk_basename",
                                          default="savpres_ss_displ.vtk")
    displVTKBaseName.meta['tip'] = "Base filename of VTK displacement output."

    displCSVFileName = pyre.inventory.str("displ_csv_filename",
                                          default="savpres_ss_displ.csv")
    displCSVFileName.meta['tip'] = "Filename for CSV displacement output."

    velVTKBaseName = pyre.inventory.str("vel_vtk_basename",
                                        default="savpres_ss_vel.vtk")
    velVTKBaseName.meta['tip'] = "Base filename of VTK velocity output."

    velCSVFileName = pyre.inventory.str("vel_csv_filename",
                                        default="savpres_ss_vel.csv")
    velCSVFileName.meta['tip'] = "Filename for CSV velocity output."

    displScaleFactor = pyre.inventory.float("displ_scale_factor", default=1.0)
    displScaleFactor.meta['tip'] = "Scale factor for displacement output."

    velScaleFactor = pyre.inventory.float("vel_scale_factor", default=1.0)
    velScaleFactor.meta['tip'] = "Scale factor for velocity output."

    coordScaleFactor = pyre.inventory.float("coord_scale_factor", default=1.0)
    coordScaleFactor.meta['tip'] = "Scale factor for output coordinates."

    coordUnits = pyre.inventory.str("coord_units", default="m")
    coordUnits.meta['tip'] = "Units used for output coordinates."

    timeUnits = pyre.inventory.dimensional("time_units", default=1.0*year)
    timeUnits.meta['tip'] = "Time units to use for output filenames."

    timeStampWidth = pyre.inventory.int("time_stamp_width", default=4)
    timeStampWidth.meta['tip'] = "Number digits in output filename time stamp."

    title = pyre.inventory.str("title",
                               default="Savage & Prescott strike-slip solution")
    title.meta['tip'] = "Title to appear at the top of VTK files."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="savpres_ss"):
    Application.__init__(self, name)
    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._genPoints()
    self._genSolution()
    self.points *= self.coordScaleFactor
    if self.outputDisplVTK:
      self._writeSolutionVTK("displacement")
    if self.outputVelVTK:
      self._writeSolutionVTK("velocity")
    if self.outputDisplCSV:
      self._writeSolutionCSV("displacement")
    if self.outputVelCSV:
      self._writeSolutionCSV("velocity")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.elasThick = self.inventory.elasThick.value
    self.lockDepth = self.inventory.lockDepth.value
    self.recurrenceTime = self.inventory.recurrenceTime.value
    self.viscosity = self.inventory.viscosity.value
    self.shearModulus = self.inventory.shearModulus.value
    self.velocity = self.inventory.plateVelocity.value/2.0
    self.numberCycles = self.inventory.numberCycles
    self.numberSteps = self.inventory.numberSteps
    self.numberTerms = self.inventory.numberTerms
    self.numberPoints = self.inventory.numberPoints
    self.deltaX = self.inventory.deltaX.value
    self.xEpsilon = self.inventory.xEpsilon.value
    self.outputDisplVTK = self.inventory.outputDisplVTK
    self.outputDisplCSV = self.inventory.outputDisplCSV
    self.outputVelVTK = self.inventory.outputVelVTK
    self.outputVelCSV = self.inventory.outputVelCSV
    self.displVTKBaseName = self.inventory.displVTKBaseName
    self.displCSVFileName = self.inventory.displCSVFileName
    self.velVTKBaseName = self.inventory.velVTKBaseName
    self.velCSVFileName = self.inventory.velCSVFileName
    self.displScaleFactor = self.inventory.displScaleFactor
    self.velScaleFactor = self.inventory.velScaleFactor
    self.coordScaleFactor = self.inventory.coordScaleFactor
    self.coordUnits = self.inventory.coordUnits
    self.timeUnits = self.inventory.timeUnits.value
    self.timeStampWidth = self.inventory.timeStampWidth
    self.title = self.inventory.title

    self.deltaT = self.recurrenceTime/self.numberSteps
    self.tauFac = 0.5*self.shearModulus/self.viscosity
    self.tau0 = self.recurrenceTime * self.tauFac

    return


  def _genPoints(self):
    """
    Create array of points for output along with series terms
    for each point.
    """
    self.points = numpy.zeros(self.numberPoints, dtype=numpy.float64)
    self.pointCoeff = numpy.zeros((self.numberPoints, self.numberTerms),
                                  dtype=numpy.float64)

    for point in range(self.numberPoints):
      self.points[point] = max(self.xEpsilon, point*self.deltaX)

      for term in range(self.numberTerms):
        n = term + 1
        self.pointCoeff[point, term] = 2.0 * self.lockDepth * \
                                       self.points[point]/ \
                                       (4.0 * n**2 * self.elasThick**2 - \
                                        self.lockDepth**2 + \
                                        self.points[point]**2)

    self.pointCoeff = numpy.arctan(self.pointCoeff)

    return

    
  def _genSolution(self):
    """
    Compute transient solution.
    """
    solutionU2 = numpy.zeros((self.numberCycles,
                              self.numberSteps + 1,
                              self.numberPoints),
                             dtype=numpy.float64)
    self.solutionUTot = numpy.zeros((self.numberCycles,
                                     self.numberSteps + 1,
                                     self.numberPoints),
                                    dtype=numpy.float64)
    solutionV2 = numpy.zeros((self.numberCycles,
                              self.numberSteps + 1,
                              self.numberPoints),
                             dtype=numpy.float64)
    self.solutionVTot = numpy.zeros((self.numberCycles,
                                     self.numberSteps + 1,
                                     self.numberPoints),
                                    dtype=numpy.float64)
    oneArray = numpy.ones(self.numberPoints, dtype=numpy.float64)

    for cycle in range(self.numberCycles):
      time = cycle * self.numberSteps * self.deltaT
      tau = time * self.tauFac
      if cycle > 0:
        solutionU2[cycle, :, :] += solutionU2[cycle - 1, :, :]
        solutionV2[cycle, :, :] += solutionV2[cycle - 1, :, :]

      for step in range(self.numberSteps + 1):
        if cycle == 0:
          solutionUT, solutionVT = self._u2A(tau)
        else:
          solutionUT, solutionVT = self._u2B(tau)

        solutionU2[cycle, step, :] += solutionUT
        solutionV2[cycle, step, :] += solutionVT
        self.solutionUTot[cycle, step, :] = solutionU2[cycle, step, :] + \
                                            time * self.velocity * oneArray
        self.solutionVTot[cycle, step, :] = self.tauFac * \
                                            solutionV2[cycle, step, :] + \
                                            self.velocity * oneArray
          
        time = time + self.deltaT
        tau = time * self.tauFac

    self.solutionUTot *= self.displScaleFactor
    self.solutionVTot *= self.velScaleFactor

    return


  def _timeCoeff(self, term, tau, aPrev, bPrev, factPrev):
    """
    Computes coefficients for term term and time tau.
    """
    if term == 0:
      factN = 1.0
      aN = 1.0 - math.exp(-tau)
      bN = (tau - aN)/self.tau0
    else:
      factN = term * factPrev
      aN = aPrev - tau**term * math.exp(-tau)/factN
      bN = bPrev - aN/self.tau0

    return aN, bN, factN
        
      
  def _u2A(self, tau):
    """
    Computes viscoelastic solution for times less than the recurrence time.
    """
    solutionU = numpy.zeros(self.numberPoints, dtype=numpy.float64)
    solutionV = numpy.zeros(self.numberPoints, dtype=numpy.float64)

    for point in range(self.numberPoints):
      solution = (-0.5 * math.pi + \
                  numpy.arctan(self.points[point]/self.lockDepth))/self.tau0
      solutionU[point] = tau * solution
      solutionV[point] = solution
      aPrev = 0.0
      bPrev = 0.0
      factPrev = 1.0
      for term in range(self.numberTerms):
        aN, bN, factN = self._timeCoeff(term, tau, aPrev, bPrev, factPrev)
        solutionU[point] -= bN * self.pointCoeff[point, term]
        solutionV[point] -= aN * self.pointCoeff[point, term]/self.tau0
        aPrev = aN
        bPrev = bN
        factPrev = factN

    solutionU *= 2.0 * self.velocity * self.recurrenceTime/math.pi
    solutionV *= 2.0 * self.velocity * self.recurrenceTime/math.pi
    return [solutionU, solutionV]
        
      
  def _u2B(self, tau):
    """
    Computes viscoelastic solution for times greater than the recurrence time.
    """
    solutionU = numpy.zeros(self.numberPoints, dtype=numpy.float64)
    solutionV = numpy.zeros(self.numberPoints, dtype=numpy.float64)
    tau2 = tau - self.tau0

    for point in range(self.numberPoints):
      a1Prev = 0.0
      b1Prev = 0.0
      fact1Prev = 1.0
      a2Prev = 0.0
      b2Prev = 0.0
      fact2Prev = 1.0
      for term in range(self.numberTerms):
        a1N, b1N, fact1N = self._timeCoeff(term, tau, a1Prev, b1Prev, fact1Prev)
        a2N, b2N, fact2N = self._timeCoeff(term, tau2, a2Prev, b2Prev,
                                           fact2Prev)
        daDt = tau2**term * math.exp(-tau2)/fact2N
        solutionU[point] += self.pointCoeff[point, term] * \
                           (b2N - b1N + a2N)
        solutionV[point] += self.pointCoeff[point, term] * \
                            (a2N/self.tau0 - a1N/self.tau0 + daDt)
        a1Prev = a1N
        b1Prev = b1N
        fact1Prev = fact1N
        a2Prev = a2N
        b2Prev = b2N
        fact2Prev = fact2N

    solutionU *= 2.0 * self.velocity * self.recurrenceTime/math.pi
    solutionV *= 2.0 * self.velocity * self.recurrenceTime/math.pi
    return [solutionU, solutionV]
    
      
  def _writeSolutionVTK(self, solutionType):
    """
    Generate VTK filename and write results to file.
    """
    
    if solutionType == "displacement":
      VTKBaseName = self.displVTKBaseName
      solution = self.solutionUTot
    else:
      VTKBaseName = self.velVTKBaseName
      solution = self.solutionVTot

    if VTKBaseName.endswith('.vtk'):
      fileBase = VTKBaseName[:VTKBaseName.rfind('.vtk')]
    elif VTKBaseName.endswith('.VTK'):
      fileBase = VTKBaseName[:VTKBaseName.rfind('.VTK')]
    else:
      fileBase = VTKBaseName

    for cycle in range(self.numberCycles):
      fileBaseCycle = fileBase + "_c" + str(cycle) + "_t"
      time = 0.0

      for step in range(self.numberSteps + 1):
        timeStampInt = int(time/self.timeUnits)
        timeStampString = repr(timeStampInt).rjust(self.timeStampWidth, '0')
        VTKFile = fileBaseCycle + timeStampString + ".vtk"
        f = open(VTKFile, 'w')
        self._writeVTK(f, solution, solutionType, cycle, step)
        f.close()
        time += self.deltaT

    return


  def _writeVTK(self, f, solution, solutionType, cycle, step):
    """
    Writes solution to VTK file as a set of points.
    """
    f.write('# vtk DataFile Version 2.0\n')
    f.write(self.title + '\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS '+str(self.numberPoints)+' double\n')
    y = 0.0
    z = 0.0
    for point in self.points:
      f.write(' %.12g   %.12g   %.12g\n' % (point, y, z))

    f.write('\n')
    f.write('POINT_DATA '+str(self.numberPoints)+'\n')
    f.write('SCALARS '+solutionType+' double 3\n')
    f.write('LOOKUP_TABLE default\n')
    uX = 0.0
    uZ = 0.0
    for point in range(self.numberPoints):
      f.write(' %.12g   %.12g   %.12g\n' %
              (uX, solution[cycle, step, point], uZ))
    
    return
    
      
  def _writeSolutionCSV(self, solutionType):
    """
    Write solution to a CSV file.
    """
    
    if solutionType == "displacement":
      CSVFileName = self.displCSVFileName
      solution = self.solutionUTot
    else:
      CSVFileName = self.velCSVFileName
      solution = self.solutionVTot

    f = open(CSVFileName, 'w')
    head = "Distance from Fault (" + self.coordUnits + ")"
    for cycle in range(self.numberCycles):
      cycleHead = "Cycle " + str(cycle) + " t = "
      time = 0.0

      for step in range(self.numberSteps + 1):
        timeStampInt = int(time/self.timeUnits)
        timeStampString = repr(timeStampInt) + " years"
        head += "," + cycleHead + timeStampString
        time += self.deltaT

    f.write('%s\n' % head)

    for point in range(self.numberPoints):
      f.write(' %.12g' % (self.points[point]))
      for cycle in range(self.numberCycles):
        for step in range(self.numberSteps + 1):
          f.write(', %.12g' % solution[cycle, step, point])
      f.write('\n')

    f.close()

    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Savpres_ss()
  app.run()

# End of file
