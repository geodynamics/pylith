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

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace bc {
    class DirichletBoundary;
    class TestDirichletBoundary; // unit testing
  } // bc

  namespace topology {
    class FieldUniform; // USES FieldUniform
    class SolutionFields; // USES SolutionFields
  } // topology
} // pylith

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
  const ALE::Obj<SieveSubMesh>& boundaryMesh(void) const;

  /** Get vertex field with BC information.
   *
   * @param fieldType Type of field.
   * @param name Name of field.
   * @param mesh Finite-element mesh.
   * @param fields Solution fields.
   *
   * @returns Field over vertices.
   */
  const topology::Field&
  vertexField(const char* name,
	      const topology::Mesh& mesh,
	      const topology::SolutionFields& fields);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Extract submesh associated with boundary.
   *
   * @param mesh Finite-element mesh.
   */
  void _createBoundaryMesh(const topology::Mesh& mesh);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  DirichletBoundary(const DirichletBoundary& m);

  /// Not implemented
  const DirichletBoundary& operator=(const DirichletBoundary& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  ALE::Obj<SieveSubMesh> _boundaryMesh; ///< Boundary mesh.
  topology::FieldUniform* _tmpField; ///< Temporary field for output.

}; // class DirichletBoundary

#include "DirichletBoundary.icc" // inline methods

#endif // pylith_bc_dirichletboundary_hh


// End of file 
