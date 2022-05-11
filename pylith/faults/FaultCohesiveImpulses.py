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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

import numpy

from .FaultCohesive import FaultCohesive
from .faults import FaultCohesiveImpulses as ModuleFaultCohesiveImpulses

def validateDOF(value):
    """Validate list of fixed degrees of freedom.
    """
    def error():
        raise ValueError("'impuluse_dof' must be a zero based list of indices of degrees of "
          "freedom at a vertex.")
    try:
        size = len(value)
        if 0 == size:
            raise error()
        num = list(map(int, value))
        for v in num:
            if v < 0:
                raise error()
    except:
        error()
    return num

class FaultCohesiveImpulses(FaultCohesive, ModuleFaultCohesiveImpulses):
    """
    Fault surface with slip impulses for Green's functions implemented with cohesive cells.

    The comopnents 

    Implements `FaultCohesiveKin`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.greensfns]
            interfaces = [fault]
            interfaces.fault = pylith.faults.FaultCohesiveImpulses

            [pylithapp.greensfns.interfaces.fault]
            label = fault
            label_value = 20
            
            # Impulses for left-lateral slip (dof=1)
            impulse_dof = [1]
            threshold = 0.5

            # Create impulses at all points on the fault by specifying a uniform amplitude of 1.0.
            # Impulses will be applied at any location with a slip component greater than the threshold.
            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Slip impulse amplitude
            db_auxiliary_field.values = [slip_left_lateral, slip_opening]
            db_auxiliary_field.data = [1.0*m, 0.0*m]

            # Represent the impulse as a linear variation in slip centered on each point.
            auxiliary_subfields.slip.basis_order = 1
            """
    }

    import pythia.pyre.inventory
    from pythia.pyre.units.length import m

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)

    threshold = pythia.pyre.inventory.dimensional("threshold", default=1.0e-6*m, validator=pythia.pyre.inventory.greaterEqual(0.0*m))
    threshold.meta['tip'] = "Threshold for non-zero amplitude."

    impulseDOF = pythia.pyre.inventory.list("impulse_dof", default=[], validator=validateDOF)
    impulseDOF.meta['tip'] = "Indices of impulse components; 0=fault opening, 1=left lateral, 2=reverse (3D only)."

    def __init__(self, name="faultcohesiveimpulses"):
        """
        Initialize configuration.
        """
        FaultCohesive.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsFault import AuxSubfieldsFault
        self.auxiliarySubfields = AuxSubfieldsFault("auxiliary_subfields")

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()

        if 0 == comm.rank:
            self._info.log(f"Pre-initializing fault '{self.labelName}={self.labelValue}'.")
        FaultCohesive.preinitialize(self, problem)
        ModuleFaultCohesiveImpulses.setThreshold(self, self.threshold.value)
        impulseDOF = numpy.array(self.impulseDOF, dtype=numpy.int32)
        ModuleFaultCohesiveImpulses.setImpulseDOF(self, impulseDOF)
  

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultCohesiveImpulses.verifyConfiguration(self, self.mesh())

    def _createModuleObj(self):
        """Create handle to C++ FaultCohesiveImpulses.
        """
        ModuleFaultCohesiveImpulses.__init__(self)
  

# FACTORIES ////////////////////////////////////////////////////////////

def fault():
  """
  Factory associated with FaultCohesiveImpulses.
  """
  return FaultCohesiveImpulses()


# End of file 
