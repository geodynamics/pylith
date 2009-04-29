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

    class FaultCohesive : public Fault
    { // class FaultCohesive

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FaultCohesive(void);
      
      /// Destructor.
      virtual
      ~FaultCohesive(void);
      
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
      
      /** Adjust mesh topology for fault implementation.
       *
       * @param mesh PETSc mesh.
       * @param flipFault Flip fault orientation.
       */
      void adjustTopology(pylith::topology::Mesh* const mesh,
			  const bool flipFault =false);
      
      // PROTECTED METHODS //////////////////////////////////////////////////
    protected :
      
      /** Cohesive cells use Lagrange multiplier constraints?
       *
       * @returns True if implementation using Lagrange multiplier
       * constraints, false otherwise.
       */
      virtual
      bool _useLagrangeConstraints(void) const = 0;
      
    }; // class FaultCohesive

  } // faults
} // pylith


// End of file 
