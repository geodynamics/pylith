#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams
#  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


from pyre.components.Component import Component


class SimulationController(Component):


    class Inventory(Component.Inventory):

        import pyre.inventory
        
        restart = pyre.inventory.bool("restart", default=False)
        restartStep = pyre.inventory.int("restartStep", default=0)
        dt = pyre.inventory.float("dt", default=1.0)
        frameFrequency = pyre.inventory.int("frameFrequency", default=100)
        checkpointFrequency = pyre.inventory.int("checkpointFrequency", default=100)
        monitoringFrequency = pyre.inventory.int("monitoringFrequency", default=100)
        outputDirectory = pyre.inventory.str("outputDirectory", default=".")
        

    def launch(self, app):
        self.solver.launch(app)
        return


    def restart(self):
        return


    def march(self, totalTime=None, steps=0):
        """explicit time loop"""

        # the main simulation time loop
        while 1:

            # notify solvers we are starting a new timestep
            self.startTimestep()

            # synchronize boundary information
            self.applyBoundaryConditions()

            # do io
            self.save()

            # compute an acceptable timestep
            dt = self.stableTimestep()

            # advance
            self.advance(dt)

            # update smulation clock and step number
            self.clock = self.clock + dt
            self.step = self.step + 1

            # notify solver we finished a timestep
            self.endTimestep()

            # are we done?
            if totalTime and self.clock >= totalTime:
                break
            if steps and self.step >= steps:
                break

        # end of time advance loop           

        # Notify solver we are done
        self.endSimulation()

        return


    def startTimestep(self):
        self.solver.newStep(self.clock, self.step)
        return

    
    def applyBoundaryConditions(self):
        self.solver.applyBoundaryConditions()
        return


    def stableTimestep(self):
        dt = self.solver.stableTimestep()
        return dt


    def advance(self, dt):
        self.solver.advance(dt)
        return


    def save(self):

        step = self.step
        directory = None
        
        if not step % self.inventory.monitoringFrequency:
            if not directory:
                directory = self._makeDirectory(step)
            self.solver.publishState(directory)
            
        # save visualization information
        if not step % self.inventory.frameFrequency:
            if not directory:
                directory = self._makeDirectory(step)
            self.solver.plotFile(directory)

        # dump checkpoint
        # don't checkpoint if we just restarted from checkpoint
        if not self.inventory.restart or step > self.inventory.restartStep:
            if step and not step % self.inventory.checkpointFrequency:
                if not directory:
                    directory = self._makeDirectory(step)
                self.solver.checkpoint(directory)
            
        return


    def endTimestep(self):
        self.solver.endTimestep(self.clock)
        return


    def endSimulation(self):
        self.solver.endSimulation(self.step, self.clock)
        return


    def _makeDirectory(self, step):
        outdirbase = self.inventory.outputDirectory + "/" + self.solver.name
        outdirformat = "%s-%%05d" % outdirbase
        directory = outdirformat % step
        try:
            import os
            os.makedirs(directory)
        except OSError, error:
            import errno
            errorCode, msg =  error

            # ignore failure if attempting to create a directory that exists
            if errorCode == errno.EEXIST:
                return directory
            
            raise error

        return directory


    def __init__(self, name=None, facility=None):
        if name is None:
            name = "controller"

        if facility is None:
            facility = "controller"
            
        Component.__init__(self, name, facility)

        from pyre.units.time import second

        self.step = 0
        self.clock = 0.0 * second
        self.solver = None

        return


# version
__id__ = "$Id: SimulationController.py,v 1.1.1.1 2005/03/08 16:13:46 aivazis Exp $"

#
# End of file
