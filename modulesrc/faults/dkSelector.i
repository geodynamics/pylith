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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/dkSelector.i 
 *
 * @brief Python interface to C++ Fault object.
 */

namespace pylith {
  namespace faults {

    class dkSelector
    { // class dkSelector

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      dkSelector(void);
      
      /// Destructor.
      ~dkSelector(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set spatial database for dkSelector
       *
       * @param db Spatial database
       */
      void dbdksel(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Initialize slip time function.
       *
       * @param faultMesh Finite-element mesh of fault.
       * @param normalizer Nondimensionalization of scales.
       */
      void initialize(const pylith::topology::SubMesh& faultMesh,
		      const spatialdata::units::Nondimensional& normalizer);
      
      /** Get dk on fault surface (time will be implemented through this guy)
       *
       * @param dk DK selector field over fault surface
       *
       * @returns dk for the time
       */
      void dk(pylith::topology::Field<pylith::topology::SubMesh>* const dkField);
  
    }; // class dkSelector

  } // faults
} // pylith


// End of file 
