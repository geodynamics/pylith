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

/** @file modulesrc/faults/TimeHistorySlipFn.i
 *
 * @brief Python interface to C++ Fault object.
 */

namespace pylith {
  namespace faults {

    class TimeHistorySlipFn : public SlipTimeFn
    { // class TimeHistorySlipFn

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      TimeHistorySlipFn(void);
      
      /// Destructor.
      ~TimeHistorySlipFn(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set spatial database for amplitude of slip time history.
       *
       * @param db Spatial database
       */
      void dbAmplitude(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set spatial database for slip initiation time.
       *
       * @param db Spatial database
       */
      void dbSlipTime(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set time history.
       *
       * @param th Time history.
       */
      void dbTimeHistory(spatialdata::spatialdb::TimeHistory* const th);

      /** Initialize slip time function.
       *
       * @param faultMesh Finite-element mesh of fault.
       * @param cs Coordinate system for mesh
       * @param normalizer Nondimensionalization of scales.
       * @param originTime Origin time for earthquake source.
       */
      void initialize(const pylith::topology::Mesh& faultMesh,
		      const spatialdata::units::Nondimensional& normalizer,
		      const PylithScalar originTime =0.0);
      
      /** Get slip on fault surface at time t.
       *
       * @param slipField Slip field over fault surface.
       * @param t Time t.
       *
       * @returns Slip vector as left-lateral/reverse/normal.
       */
      void slip(pylith::topology::Field* const slipField,
		const PylithScalar t);
  
      /** Get final slip.
       *
       * @returns Final slip.
       */
      const pylith::topology::Field& finalSlip(void);
      
      /** Get time when slip begins at each point.
       *
       * @returns Time when slip begins.
       */
      const pylith::topology::Field& slipTime(void);

    }; // class TimeHistorySlipFn

  } // faults
} // pylith


// End of file 
