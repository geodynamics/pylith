#!/usr/bin/env nemesis
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

__requires__ = "PyLith"


# ======================================================================
from pyre.applications.Script import Script
class ParametersApp(Script):
  """
  Application for printing current PyLith parameters to a text file.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Script.Inventory):
    """
    Python object for managing ParametersApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ParametersApp facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Filename for output.
    ## @li \b verbose Set to true to print location where each
    ##   parameter is set.
    ##
    ## \b Facilities
    ## @li none

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="parameters.txt")
    filename.meta['tip'] = "Filename for output."

    verbose = pyre.inventory.bool("verbose", default=False)
    verbose.meta['tip'] = "Set to true to print location where " \
        "each parameter is set."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pylithinfo"):
    """
    Constructor.
    """
    Script.__init__(self, name)
    self._tab = "  "
    self._printPropertyFn = self._printPropertyBasic
    self._printFacilityFn = self._printFacilityBasic
    return


  def main(self, *args, **kwds):
    """
    Main entry point for application.
    """
    
    from pylith.apps.PyLithApp import InfoApp
    targetapp = InfoApp()
    targetapp.run()

    if self.verbose:
      self._printPropertyFn = self._printPropertyVerbose
      self._printFacilityFn = self._printFacilityVerbose
    else:
      self._printPropertyFn = self._printPropertyBasic
      self._printFacilityFn = self._printFacilityBasic

    depth = 0
    fout = open(self.filename, "w")
    fout.write("Application: %s %s\n" % (targetapp.name, targetapp))
    self._printParams(fout, targetapp, depth+1)
    fout.close()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members using inventory.
    """
    Script._configure(self)
    self.filename = self.inventory.filename
    self.verbose = self.inventory.verbose
    return


  def _printParams(self, fout, obj, depth):
    """
    Print objects parameters to fout.
    """
    propertyNames = obj.inventory.propertyNames()
    facilityNames = obj.inventory.facilityNames()

    propertiesOmit = ["help", 
                      "help-components",
                      "help-persistence",
                      "help-properties",
                      "typos",
                      ]
    facilitiesOmit = ["weaver",
                      ]

    propertyNames.sort()
    facilityNames.sort()

    for name in propertyNames:
      if name in facilityNames or name in propertiesOmit:
        continue
      trait = obj.inventory.getTrait(name)
      descriptor = obj.inventory.getTraitDescriptor(name)
      self._printPropertyFn(fout, name, trait, descriptor, depth)
    for name in facilityNames:
      if name in facilitiesOmit:
        continue
      trait = obj.inventory.getTrait(name)
      descriptor = obj.inventory.getTraitDescriptor(name)
      self._printFacilityFn(fout, name, trait, descriptor, depth)
    return


  def _printPropertyBasic(self, fout, name, trait, descriptor, depth):
    """
    Print property name, type, and value.
    """
    indent = self._tab*depth
    fout.write("%s%s (%s) = %s\n" % \
                 (indent, name, trait.type, descriptor.value))
    return


  def _printFacilityBasic(self, fout, name, trait, descriptor, depth):
    """
    Print facility name, type, and value.
    """
    indent = self._tab*depth
    fout.write("%s%s = %s (%s)\n" % \
                 (indent, name, descriptor.value.name, descriptor.value))
    self._printParams(fout, descriptor.value, depth+1)
    return


  def _printPropertyVerbose(self, fout, name, trait, descriptor, depth):
    """
    Print property, name, type, value, description, and location set.
    """
    indent = self._tab*depth
    fout.write("\n%s%s (%s) = %s\n" % \
                 (indent, name, trait.type, descriptor.value))    

    indent += self._tab
    try:
      description = trait.meta['tip']
    except KeyError:
      description = "No description available."
    fout.write("%sDescription: %s\n" % (indent, description))
    fout.write("%sSet from: %s\n" % (indent, descriptor.locator))
    return


  def _printFacilityVerbose(self, fout, name, trait, descriptor, depth):
    """
    Print facility name, type, value, description, and location set.
    """
    indent = self._tab*depth
    fout.write("\n%s%s = %s (%s)\n" % \
                 (indent, name, descriptor.value.name, descriptor.value))

    indent += self._tab
    try:
      description = trait.meta['tip']
    except KeyError:
      description = "No description available."
    fout.write("%sDescription: %s\n" % (indent, description))
    fout.write("%sSet from: %s\n" % (indent, descriptor.locator))
    fout.write("%sConfigurable as: %s\n" % \
                 (indent, ", ".join(descriptor.value.aliases)))
    
    self._printParams(fout, descriptor.value, depth+1)
    return


# ----------------------------------------------------------------------
if __name__ == "__main__":

  app = ParametersApp()
  app.run()


# End of file 
