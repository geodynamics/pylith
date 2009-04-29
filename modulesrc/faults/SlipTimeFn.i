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
		      const double originTime =0.0) = 0;

      /** Get slip on fault surface at time t.
       *
       * @param slipField Slip field over fault surface.
       * @param t Time t.
       *
       * @returns Slip vector as left-lateral/reverse/normal.
       */
      virtual
      void slip(pylith::topology::Field<pylith::topology::SubMesh>* const slipField,
		const double t) = 0;
  
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
		    const double t0,
		    const double t1) = 0;
      
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

    }; // class SlipTimeFn

  } // faults
} // pylith


// End of file 
