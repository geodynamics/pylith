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

/** @file modulesrc/bc/BoundaryCondition.i
 *
 * @brief Python interface to C++ BoundaryCondition object.
 */

namespace pylith {
  namespace bc {

    class BoundaryCondition
    { // class BoundaryCondition

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      BoundaryCondition(void);

      /// Destructor.
      virtual
      ~BoundaryCondition(void);

      /** Set label of boundary condition surface.
       *
       * @param value Label of surface (from mesh generator).
       */
      void label(const char* value);

      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      const char* label(void) const;

      /** Set database for boundary condition parameters.
       *
       * @param db Spatial database
       */
      void db(spatialdata::spatialdb::SpatialDB* const db);

      /** Verify configuration.
       *
       * @param mesh Finite-element mesh.
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

      /** Initialize boundary condition.
       *
       * @param mesh Finite-element mesh.
       * @param upDir Vertical direction (somtimes used in 3-D problems).
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      const double upDir[3]) = 0;

    }; // class BoundaryCondition

  } // bc
} // pylith


// End of file 
