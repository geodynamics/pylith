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

/** @file modulesrc/faults/SlipTimeFn.i
 *
 * @brief Python interface to C++ Fault object.
 */

namespace pylith {
  namespace faults {

    class SlipTimeFn
    { // class SlipTimeFn

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      SlipTimeFn(void);

      /// Destructor.
      virtual
      ~SlipTimeFn(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Initialize slip time function.
       *
       * @param faultMesh Finite-element mesh of fault.
       * @param cs Coordinate system for mesh
       * @param normalizer Nondimensionalization of scales.
       * @param originTime Origin time for earthquake source.
       */
      virtual
      void initialize(const pylith::topology::SubMesh& faultMesh,
		      const spatialdata::units::Nondimensional& normalizer,
		      const PylithScalar originTime =0.0) = 0;

      /** Get slip on fault surface at time t.
       *
       * @param slipField Slip field over fault surface.
       * @param t Time t.
       *
       * @returns Slip vector as left-lateral/reverse/normal.
       */
      virtual
      void slip(pylith::topology::Field<pylith::topology::SubMesh>* const slipField,
		const PylithScalar t) = 0;
  
      /** Get slip increment on fault surface between time t0 and t1.
       *
       * @param slipField Slip field over fault surface.
       * @param t0 Time t.
       * @param t1 Time t+dt.
       * 
       * @returns Increment in slip vector as left-lateral/reverse/normal.
       */
      virtual
      void slipIncr(pylith::topology::Field<pylith::topology::SubMesh>* slipField,
		    const PylithScalar t0,
		    const PylithScalar t1) = 0;
      
      /** Get final slip.
       *
       * @returns Final slip.
       */
      virtual
      const pylith::topology::Field<pylith::topology::SubMesh>& finalSlip(void) = 0;

      /** Get time when slip begins at each point.
       *
       * @returns Time when slip begins.
       */
      virtual
      const pylith::topology::Field<pylith::topology::SubMesh>& slipTime(void) = 0;

      /** Get parameter fields.
       *
       * @returns Parameter fields.
       */
      const pylith::topology::Fields<pylith::topology::Field<pylith::topology::SubMesh> >*
      parameterFields(void) const;

    }; // class SlipTimeFn

  } // faults
} // pylith


// End of file 
