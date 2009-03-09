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

/** @file modulesrc/feassemble/Constraint.i
 *
 * @brief Python interface to C++ abstract base Constraint.
 */

namespace pylith {
  namespace feassemble {

    class Constraint
    { // class Constraint

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      Constraint(void);

      /// Destructor.
      virtual
      ~Constraint(void);

      /** Set manager of scales used to nondimensionalize problem.
       *
       * @param dim Nondimensionalizer.
       */
      void normalizer(const spatialdata::units::Nondimensional& dim);

      /** Set number of degrees of freedom that are constrained at
       * points in field.
       *
       * @param field Solution field
       */
      virtual
      void setConstraintSizes(const pylith::topology::Field<pylith::topology::Mesh>& field) = 0;

      /** Set which degrees of freedom are constrained at points in field.
       *
       * @param field Solution field
       */
      virtual
      void setConstraints(const pylith::topology::Field<pylith::topology::Mesh>& field) = 0;

      /** Set flag for setting constraints for total field solution or
       *  incremental field solution.
       *
       * @param flag True if using incremental solution, false otherwise.
       */
      virtual
      void useSolnIncr(const bool flag);
      
      /** Set values in field.
       *
       * @param t Current time
       * @param field Solution field
       */
      virtual
      void setField(const double t,
		    const pylith::topology::Field<pylith::topology::Mesh>& field) = 0;
      
    }; // class Constraint

  } // feassemble
} // pylith


// End of file 
