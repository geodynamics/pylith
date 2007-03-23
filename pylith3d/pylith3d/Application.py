#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from cig.cs.petsc import PetscApplication


class Application(PetscApplication):


    name = "pylith3d"


    # Tell the framework where to find PETSc functions.
    import pylith3d as petsc


    # Use PETSc-style command line parsing.
    from cig.cs.petsc import PetscCommandlineParser as CommandlineParser


    def readSamplePoints(self, filename):
        '''Read in the sampling locations
        - One point per line, three values per line (x,y,z)
        - Returns a Numeric array'''
        import Numeric
        f = file(filename)
        points = []
        for line in f.readlines():
            points.append([float(v) for v in line.strip().split(' ')])
        f.close()
        return Numeric.array(points)


    def outputSampleValues(self, filename, values):
        '''sample# sample values impluse# impulse type'''
        # Computing normal to the fault:
        #   Split nodes define the fault
        #   Get all fault faces for a node
        #   Area weighted average of normals
        f = file(filename, 'w')
        for v, values in enumerate(values):
            write(f, '%d %g %g %g 1 0' % (v, values[0], values[1], values[2]))
        f.close()
        return


    def greenFunction(self, points):
        """
        # Beginning of loop that loops over split node sets, creating
        # an 'impulse' for each one and outputting response values.
        # Below at present is a quasi-C version of the needed code.
SectionReal splitField;

# Need bindings for this
ierr = MeshGetSectionPair(mesh, "split", &splitField);
// Loop over split nodes
for() {
  // Loop over elements
  for() {
# Need bindings for this
    ierr = SectionPairSetFiberDimension(splitField, e, 1);
  }
# Need bindings for this
  ierr = SectionPairAllocate(splitField);
  // Loop over elements
  for() {
    PetscPair value;

    value.i = node;
    value.x = ;
    value.y = ;
    value.z = ;
# Need bindings for this
    ierr = SectionPairUpdate(splitField, e, &value);
# Major problem right now:  This just updates PETSc/Sieve's copy of splitField.
# It does not change the values within PyLith, which have been read from
# per-process input files.
  }
  // Solve
        pl3drun.solveElastic()
# Need bindings for this
  ierr = SectionPairClear(splitField);
}
"""
        scanner = self.inventory.scanner
        values = scanner.interpolatePoints(points)
        self.outputSampleValues(scanner.fileRoot+'.output', values)
        return


    def main(self, *args, **kwds):
    
#        from time import clock as now
#        start = now()
        import pylith3d

        scanner = self.inventory.scanner
        green = self.inventory.green
        
        if green:
            points      = self.readSamplePoints(scanner.macroString(scanner.metainventory.sampleLocationFile))

        mesh = pylith3d.processMesh(scanner.macroString(scanner.metainventory.bcInputFile),
                                    scanner.macroString(scanner.metainventory.inputFileRoot),
                                    scanner.inventory.interpolateMesh,
                                    scanner.inventory.partitioner)

        scanner.mesh = mesh

        scanner.initialize()
        
        scanner.setup()
        scanner.read()
        scanner.numberequations()
        scanner.sortmesh()
        scanner.sparsesetup()
        scanner.allocateremaining()
        scanner.meshwrite()

        if green:
            self.greenFunction(points)
        else:
            scanner.run()
#        finish = now()
#        usertime = finish - start
#        print "Total user time:  %g" % usertime
        return


    class Inventory(PetscApplication.Inventory):

        import pyre.inventory
        from cig.cs.petsc import PetscProperty
        from Pylith3d_scan import Pylith3d_scan

        green = pyre.inventory.bool("green")

        scanner = pyre.inventory.facility("scanner", factory=Pylith3d_scan)

        # declare PETSc options that are of interest to PyLith
        ksp_monitor        = PetscProperty()
        ksp_view           = PetscProperty()
        ksp_rtol           = PetscProperty()
        log_summary        = PetscProperty()
        pc_type            = PetscProperty()
        sub_pc_type        = PetscProperty()
        start_in_debugger  = PetscProperty()
        debugger_pause     = PetscProperty()


# version
# $Id: Application.py,v 1.5 2005/04/15 00:18:21 willic3 Exp $

# End of file 
