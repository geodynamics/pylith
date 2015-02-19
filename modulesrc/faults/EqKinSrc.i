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
      void originTime(const PylithScalar value);
      
      /** Get origin time for earthquake source.
       *
       * @returns Origin time for earthquake source.
       */
      PylithScalar originTime(void) const;
      
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
      void initialize(const pylith::topology::Mesh& faultMesh,
		      const spatialdata::units::Nondimensional& normalizer);

      /** Get slip on fault surface at time t.
       *
       * @param slipField Slip field over fault mesh.
       * @param t Time t.
       */
      void slip(pylith::topology::Field* const slipField,
		const PylithScalar t);

      /** Get final slip.
       *
       * @returns Final slip.
       */
      const pylith::topology::Field& finalSlip(void) const;
      
      /** Get time when slip begins at each point.
       *
       * @returns Time when slip begins.
       */
      const pylith::topology::Field& slipTime(void) const;

    }; // class EqKinSrc

  } // faults
} // pylith


// End of file 
