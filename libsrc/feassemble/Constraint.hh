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

/** @file libsrc/feassemble/Constraint.hh
 *
 * @brief C++ abstract base class defining interface for constraints
 * applied to finite-elements.
 */

#if !defined(pylith_feassemble_constraint_hh)
#define pylith_feassemble_constraint_hh

// Include directives ---------------------------------------------------

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace feassemble {
    class Constraint;
    class TestConstraint; // unit testing
  } // feassemble

  namespace topology {
    class Field; // USES Field'
    class Mesh; // USES Mesh
  } // feassemble
} // pylith

namespace spatialdata {
  namespace units {
    class Nondimensional; // USES Nondimensional
  } // units
} // spatialdata

// Constraint -----------------------------------------------------------
class pylith::feassemble::Constraint
{ // class Constraint
  friend class TestConstraint; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
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

  /** Set number of degrees of freedom that are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh Finite-element mesh.
   */
  virtual
  void setConstraintSizes(const topology::Field& field,
			  const topology::Mesh& mesh) = 0;

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh Finite-element mesh.
   */
  virtual
  void setConstraints(const topology::Field& field,
		      const topology::Mesh& mesh) = 0;

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
   * @param mesh Finite-element mesh.
   */
  virtual
  void setField(const double t,
		const topology::Field& field,
		const topology::Mesh& mesh) = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer.

  /// Flag indicating whether to set constraints for a total field
  /// solution or an incremental field solution
  bool _useSolnIncr;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Constraint(const Constraint& m);

  /// Not implemented
  const Constraint& operator=(const Constraint& m);

}; // class Constraint

#include "Constraint.icc" // inline methods

#endif // pylith_feassemble_constraint_hh


// End of file 
