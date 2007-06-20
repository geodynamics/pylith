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

#include "pylith/utils/sievetypes.hh" // USES real_section_type

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class Constraint;
    class TestConstraint; // unit testing
  } // feassemble
} // pylith


/// C++ abstract base class defining interface for constraints applied
/// to finite-elements.
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

  /** Set number of degrees of freedom that are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  virtual
  void setConstraintSizes(const ALE::Obj<real_section_type>& field,
			  const ALE::Obj<ALE::Mesh>& mesh) = 0;

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  virtual
  void setConstraints(const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh) = 0;

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
   * @param mesh PETSc mesh
   */
  virtual
  void setField(const double t,
		const ALE::Obj<real_section_type>& field,
		const ALE::Obj<ALE::Mesh>& mesh) = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

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
