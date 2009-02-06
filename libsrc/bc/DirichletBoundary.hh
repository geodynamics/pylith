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

/** @file libsrc/bc/DirichletBoundary.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions with points on a boundary.
 */

#if !defined(pylith_bc_dirichletboundary_hh)
#define pylith_bc_dirichletboundary_hh

// Include directives ---------------------------------------------------
#include "DirichletBC.hh" // ISA DirichletBC

// DirichletBoundary ----------------------------------------------------
class pylith::bc::DirichletBoundary : public DirichletBC
{ // class DirichletBoundary
  friend class TestDirichletBoundary; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  DirichletBoundary(void);

  /// Destructor.
  ~DirichletBoundary(void);

  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Vertical direction (somtimes used in 3-D problems).
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3]);

  /** Get boundary mesh.
   *
   * @return Boundary mesh.
   */
  const topology::SubMesh& boundaryMesh(void) const;

  /** Get vertex field with BC information.
   *
   * @param name Name of field.
   * @param fields Solution fields.
   *
   * @returns Field over vertices.
   */
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields& fields);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  DirichletBoundary(const DirichletBoundary& m);

  /// Not implemented
  const DirichletBoundary& operator=(const DirichletBoundary& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  topology::SubMesh* _boundaryMesh; ///< Boundary mesh.
  topology::Field<topology::SubMesh>* _tmpField; ///< Temporary field for output.

}; // class DirichletBoundary

#include "DirichletBoundary.icc" // inline methods

#endif // pylith_bc_dirichletboundary_hh


// End of file 
