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

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Get parameter fields.
       *
       * @returns Parameter fields.
       */
      const pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >*
      parameterFields(void) const;

    }; // class BoundaryConditionPoints

  } // bc
} // pylith


// End of file 
