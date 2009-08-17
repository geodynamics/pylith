// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
			  public pylith::feassemble::Integrator<pylith::feassemble::Quadrature<pylith::topology::SubMesh> >
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
       * @param flag True if using fault mesh, false if using vertices.
       */
      void useFaultMesh(const bool flag);
      
      // TEMPORARY
      /** Set filename of UCD file for fault mesh.
       *
       * @param filename Filename for UCD file.
       */
      void faultMeshFilename(const char* filename);
      
      /** Get the number of vertices on the fault.
       *
       * @param mesh PETSc mesh
       * @return Number of vertices on the fault.
       */
      int numVertices(const topology::Mesh& mesh) const;

      /** Adjust mesh topology for fault implementation.
       *
       * @param mesh PETSc mesh.
       * @param flipFault Flip fault orientation.
       */
      %apply int *INOUT {int *firstFaultVertex, int *firstFaultCell};
      void adjustTopology(pylith::topology::Mesh* const mesh,
                          int *firstFaultVertex,
                          int *firstFaultCell,
                          const bool flipFault = false);
      %clear int *firstFaultVertex, int *firstFaultCell;
      
      /** Cohesive cells use Lagrange multiplier constraints?
       *
       * @returns True if implementation using Lagrange multiplier
       * constraints, false otherwise.
       */
      virtual
      bool useLagrangeConstraints(void) const = 0;
      
    }; // class FaultCohesive

  } // faults
} // pylith


// End of file 
