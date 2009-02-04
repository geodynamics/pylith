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
    ## @li \b vtk_basename Base name for VTK output files.
    ## @li \b csv_filename Filename for CSV output.
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

    VTKBaseName = pyre.inventory.str("vtk_basename", default="savpres_ss.vtk")
    VTKBaseName.meta['tip'] = "Base filename of VTK output."

    CSVFileName = pyre.inventory.str("csv_filename", default="savpres_ss.csv")
    CSVFileName.meta['tip'] = "Filename for CSV output."

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
    self._writeSolutionVTK()
    self._writeSolutionCSV()
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
    self.VTKBaseName = self.inventory.VTKBaseName
    self.CSVFileName = self.inventory.CSVFileName
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
    self.points = numpy.zeros(self.numberPoints, dtype=float)
    self.pointCoeff = numpy.zeros((self.numberPoints, self.numberTerms),
                                  dtype=float)

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
                             dtype=float)
    self.solutionUTot = numpy.zeros((self.numberCycles,
                                     self.numberSteps + 1,
                                     self.numberPoints),
                                    dtype=float)
    oneArray = numpy.ones(self.numberPoints, dtype=float)

    for cycle in range(self.numberCycles):
      time = cycle * self.numberSteps * self.deltaT
      tau = time * self.tauFac
      if cycle > 0:
        solutionU2[cycle, :, :] += solutionU2[cycle - 1, :, :]

      for step in range(self.numberSteps + 1):
        if cycle == 0:
          solutionT = self._u2A(tau)
        else:
          solutionT = self._u2B(tau)

        solutionU2[cycle, step, :] += solutionT
        self.solutionUTot[cycle, step, :] = solutionU2[cycle, step, :] + \
                                            time * self.velocity * oneArray
          
        time = time + self.deltaT
        tau = time * self.tauFac

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
    solution = numpy.zeros(self.numberPoints, dtype=float)

    for point in range(self.numberPoints):
      solution[point] = tau*(-0.5 * math.pi + \
                             numpy.arctan(self.points[point]/self.lockDepth))/ \
                             self.tau0
      aPrev = 0.0
      bPrev = 0.0
      factPrev = 1.0
      for term in range(self.numberTerms):
        aN, bN, factN = self._timeCoeff(term, tau, aPrev, bPrev, factPrev)
        solution[point] -= bN * self.pointCoeff[point, term]
        aPrev = aN
        bPrev = bN
        factPrev = factN

    solution *= 2.0 * self.velocity * self.recurrenceTime/math.pi
    return solution
        
      
  def _u2B(self, tau):
    """
    Computes viscoelastic solution for times greater than the recurrence time.
    """
    solution = numpy.zeros(self.numberPoints, dtype=float)
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
        solution[point] += self.pointCoeff[point, term] * \
                           (b2N - b1N + a2N)
        a1Prev = a1N
        b1Prev = b1N
        fact1Prev = fact1N
        a2Prev = a2N
        b2Prev = b2N
        fact2Prev = fact2N

    solution *= 2.0 * self.velocity * self.recurrenceTime/math.pi
    return solution
    
      
  def _writeSolutionVTK(self):
    """
    Generate VTK filename and write results to file.
    """
    
    if self.VTKBaseName.endswith('.vtk'):
      fileBase = self.VTKBaseName[:self.VTKBaseName.rfind('.vtk')]
    elif self.VTKBaseName.endswith('.VTK'):
      fileBase = self.VTKBaseName[:self.VTKBaseName.rfind('.VTK')]
    else:
      fileBase = self.VTKBaseName

    for cycle in range(self.numberCycles):
      fileBaseCycle = fileBase + "_c" + str(cycle) + "_t"
      time = 0.0

      for step in range(self.numberSteps + 1):
        timeStampInt = int(time/self.timeUnits)
        timeStampString = repr(timeStampInt).rjust(self.timeStampWidth, '0')
        VTKFile = fileBaseCycle + timeStampString + ".vtk"
        f = open(VTKFile, 'w')
        self._writeVTK(f, cycle, step)
        f.close()
        time += self.deltaT

    return


  def _writeVTK(self, f, cycle, step):
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
    f.write('SCALARS displacement double 3\n')
    f.write('LOOKUP_TABLE default\n')
    uX = 0.0
    uZ = 0.0
    for point in range(self.numberPoints):
      f.write(' %.12g   %.12g   %.12g\n' %
              (uX, self.solutionUTot[cycle, step, point], uZ))
    
    return
    
      
  def _writeSolutionCSV(self):
    """
    Write solution to a CSV file.
    """
    
    f = open(self.CSVFileName, 'w')
    head = "X,Y,Z"
    for cycle in range(self.numberCycles):
      cycleHead = "c" + str(cycle) + "_t"
      time = 0.0

      for step in range(self.numberSteps + 1):
        timeStampInt = int(time/self.timeUnits)
        timeStampString = repr(timeStampInt).rjust(self.timeStampWidth, '0')
        head += "," + cycleHead + timeStampString
        time += self.deltaT

    f.write('%s\n' % head)
    y = 0.0
    z = 0.0

    for point in range(self.numberPoints):
      f.write(' %.12g, %.12g, %.12g' % (self.points[point], y, z))
      for cycle in range(self.numberCycles):
        for step in range(self.numberSteps + 1):
          f.write(', %.12g' % self.solutionUTot[cycle, step, point])
      f.write('\n')

    f.close()

    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Savpres_ss()
  app.run()

# End of file
