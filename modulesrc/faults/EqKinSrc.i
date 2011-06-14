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

/** @file modulesrc/faults/EqKinSrc.i
 *
 * @brief Python interface to C++ EqKinSrc object.
 */

namespace pylith {
  namespace faults {

    class EqKinSrc
    { // class EqKinSrc

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      EqKinSrc(void);
      
      /// Destructor.
      ~EqKinSrc(void);
      
      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set origin time for earthquake source.
       *
       * @param value Origin time for earthquake source.
       */
      void originTime(const double value);
      
      /** Get origin time for earthquake source.
       *
       * @returns Origin time for earthquake source.
       */
      double originTime(void) const;
      
      /** Set slip time function.
       *
       * @param slipfn Slip time function.
       */
      void slipfn(SlipTimeFn* slipfn);
      
      /** Initialize slip time function.
       *
       * @param faultMesh Finite-element mesh of fault.
       * @param normalizer Nondimensionalization of scales.
       */
      void initialize(const pylith::topology::SubMesh& faultMesh,
		      const spatialdata::units::Nondimensional& normalizer);

      /** Get slip on fault surface at time t.
       *
       * @param slipField Slip field over fault mesh.
       * @param t Time t.
       */
      void slip(pylith::topology::Field<pylith::topology::SubMesh>* const slipField,
		const double t);

      /** Get increment of slip on fault surface between time t0 and t1.
       *
       * @param slipField Slip increment field over fault mesh.
       * @param t0 Time for start of slip increment.
       * @param t1 Time for end of slip increment.
       */
      void slipIncr(pylith::topology::Field<pylith::topology::SubMesh>* const slipField,
		    const double t0,
		    const double t1);
      
      /** Get final slip.
       *
       * @returns Final slip.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>& finalSlip(void) const;
      
      /** Get time when slip begins at each point.
       *
       * @returns Time when slip begins.
       */
      const pylith::topology::Field<pylith::topology::SubMesh>& slipTime(void) const;

    }; // class EqKinSrc

  } // faults
} // pylith


// End of file 
