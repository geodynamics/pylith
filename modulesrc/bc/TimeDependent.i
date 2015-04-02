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

/** @file modulesrc/bc/TimeDependent.i
 *
 * @brief Python interface to C++ TimeDependent object.
 */

namespace pylith {
  namespace bc {

    class TimeDependent
    { // class TimeDependent

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      TimeDependent(void);

      /// Destructor.
      ~TimeDependent(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set database for initial values.
       *
       * @param db Spatial database
       */
      void dbInitial(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set database for rate of change of values.
       *
       * @param db Spatial database
       */
      void dbRate(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set database for change in values.
       *
       * @param db Spatial database
       */
      void dbChange(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set database for temporal evolution of change in value.
       *
       * @param db Time history database.
       */
      void dbTimeHistory(spatialdata::spatialdb::TimeHistory* const db);
      
      /** Verify configuration is acceptable.
       *
       * @param mesh Finite-element mesh
       */
      virtual
      void verifyConfiguration(const pylith::topology::Mesh& mesh) const;

      // PROTECTED METHODS //////////////////////////////////////////////////
    protected :
      
      /** Get label of boundary condition surface.
       *
       * @returns Label of surface (from mesh generator).
       */
      virtual
      const char* _getLabel(void) const = 0;
      
    }; // class TimeDependent

  } // bc
} // pylith


// End of file 
