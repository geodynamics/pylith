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

/** @file modulesrc/bc/BoundaryConditionPoints.i
 *
 * @brief Python interface to C++ BoundaryConditionPoints object.
 */

namespace pylith {
  namespace bc {

    class pylith::bc::BoundaryConditionPoints : public BoundaryCondition
    { // class BoundaryCondition

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      BoundaryConditionPoints(void);
      
      /// Destructor.
      virtual
      ~BoundaryConditionPoints(void);

    }; // class BoundaryConditionPoints

  } // bc
} // pylith


// End of file 
