// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/FaultCohesive.i
 *
 * @brief Python interface to C++ FaultCohesive object.
 */

namespace pylith {
  namespace faults {

    class FaultCohesive : public Fault,
			  public pylith::feassemble::Integrator
    { // class FaultCohesive

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesive(void);
      
      /// Destructor.
      virtual
      ~FaultCohesive(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set flag for using fault mesh or group of vertices to define
       * fault surface.
       *
       * This method is part of a KLUDGE to allow creation of cohesive
       * cells in cases where domain cells have more than one face
       * (edge for 2-D problems) on the fault.
       *
       * @param flag True if using fault mesh, false if using vertices.
       */
      void useFaultMesh(const bool flag);
      
      /** Get the number of vertices associated with the fault (before
       * fault mesh exists).
       *
       * @param mesh PETSc mesh
       * @return Number of vertices on the fault.
       */
      int numVerticesNoMesh(const pylith::topology::Mesh& mesh) const;

      /** Adjust mesh topology for fault implementation.
       *
       * @param mesh PETSc mesh.
       */
      %apply int *INOUT {int *firstFaultVertex, int *firstLagrangeVertex, int *firstFaultCell};
      void adjustTopology(pylith::topology::Mesh* const mesh,
                          int *firstFaultVertex,
                          int *firstLagrangeVertex,
                          int *firstFaultCell);
      %clear int *firstFaultVertex, int *firstLagrangeVertex, int *firstFaultCell;
      
      /** Cohesive cells use Lagrange multiplier constraints?
       *
       * @returns True if implementation using Lagrange multiplier
       * constraints, false otherwise.
       */
      bool useLagrangeConstraints(void) const;
      
      /** Get fields associated with fault.
       *
       * @returns Fields associated with fault.
       */
      const pylith::topology::Fields* fields(void) const;

    }; // class FaultCohesive

  } // faults
} // pylith


// End of file 
