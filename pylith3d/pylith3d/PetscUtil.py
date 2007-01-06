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


from pyre.applications.CommandlineParser import CommandlineParser


class PetscCommandlineParser(CommandlineParser):
    """A parser which mimics PETSc's command line processing."""

    # The logic used here is derived from the 'while' loop at the end
    # of PetscOptionsInsert().  However, this doesn't check for "bad"
    # MPICH options, as these should have been removed by MPI_Init().


    def _parse(self, argv, root):
        
        self.action = None
        self.argv = argv
        self.processed = []
        self.unprocessed = []
        
        while self.argv:
            
            arg = self.argv.pop(0)
            
            iname = self._filterNonOptionArgument(arg)
            if iname is None:
                continue
            
            if iname.lower() == "options_file":
                # NYI
                if self.argv:
                    filename = self.argv.pop(0)
                else:
                    pass # error
                continue

            if (not self.argv) or self._isOptionArgument(self.argv[0]):
                iname, value = self._parseArgument(iname)
            else:
                value = self.argv.pop(0)

            self._processArgument(iname, value, root)

        return


    def _optionPrefix(self, arg):
        for prefix in self.prefixes:
            if arg.startswith(prefix):
                return prefix
        return None


    def _isOptionArgument(self, arg):
        import string
        prefix = self._optionPrefix(arg)
        if prefix is not None:
            candidate = arg[len(prefix):]
            if (prefix == "-" and
                len(candidate) > 0 and
                candidate[0] in string.digits):
                return False
            return True
        return False


    def _filterNonOptionArgument(self, arg):
        
        prefix = self._optionPrefix(arg)
        
        if prefix is not None:
            self._debug.line("    prefix: '%s starts with '%s'" % (arg, prefix))
            candidate = arg[len(prefix):]
            return candidate
        
        # prefix matching failed; leave this argument alone
        self._debug.line("    prefix: '%s' is not an option" % arg)
        self.processed.append(arg)
        return None



from pyre.components import Component


class Petsc(Component):


    def updateConfiguration(self, registry):
        self.options = [
            (name, descriptor.value) for name, descriptor in registry.properties.iteritems()
            ]
        return []


    def getArgs(self):
        args = []
        for iname, value in self.options:
            args.append('-' + iname)
            if value != 'true':
                args.append(value)
        return args


    def __init__(self, name):
        Component.__init__(self, name, name)
        self.options = []
        return



from pyre.inventory.Facility import Facility


class PetscFacility(Facility):


    def __init__(self, name):
        Facility.__init__(self, name=name, factory=Petsc, args=[name])
        return


    def _retrieveComponent(self, instance, componentName):
        petsc = Petsc(componentName)

        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('built-in')

        return petsc, locator



# End of file 
