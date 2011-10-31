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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/faults/LiuCosSlipFn.i
 *
 * @brief Python interface to C++ Fault object.
 */

namespace pylith {
  namespace faults {

    class LiuCosSlipFn : public SlipTimeFn
    { // class LiuCosSlipFn

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      LiuCosSlipFn(void);
      
      /// Destructor.
      ~LiuCosSlipFn(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set spatial database for final slip.
       *
       * @param db Spatial database
       */
      void dbFinalSlip(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set spatial database for slip initiation time.
       *
       * @param db Spatial database
       */
      void dbSlipTime(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Set spatial database for rise time. The rise time is the time it
       * takes for the slip to increase from 0.0 to 0.95 of the final
       * value.
       *
       * @param db Spatial database
       */
      void dbRiseTime(spatialdata::spatialdb::SpatialDB* const db);
      
      /** Initialize slip time function.
       *
       * @param faultMesh Finite-element mesh of fault.
       * @param cs Coordinate system for mesh
       * @param normalizer Nondimensionalization of scales.
       * @param originTime Origin time for earthquake source.
       */
      void initialize(const pylith::topology::SubMesh& faultMesh,
		      const spatialdata::units::Nondimensional& normalizer,
		      const PylithScalar originTime =0.0);
      
      /** Get slip on fault surface at time t.
       *
       * @param slipField Slip field over fault surface.
       * @param t Time t.
       *
       * @returns Slip vector as left-lateral/reverse/normal.
       */
      void slip(pylith::topology::Field<pylith::topology::SubMesh>* const slipField,
		const PylithScalar t);
  
      /** Get slip increment on fault surface between time t0 and t1.
       *
       * @param slipField Slip field over fault surface.
       * @param t0 Time t.
       * @param t1 Time t+dt.
       * 
       * @returns Increment in slip vector as left-lateral/reverse/normal.
       */
      void slipIncr(pylith::topology::Field<pylith::topology::SubMesh>* slipField,
		    const PylithScalar t0,
		    const PylithScalar t1);
      
      /** Get final slip.
       *
       * @returns Final slip.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>& finalSlip(void);
      
      /** Get time when slip begins at each point.
       *
       * @returns Time when slip begins.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>& slipTime(void);

    }; // class LiuCosSlipFn

  } // faults
} // pylith


// End of file 
