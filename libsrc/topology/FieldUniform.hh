// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file pylith/topology/FieldUniform.hh
 *
 * @brief Vector field with uniform fiber dimension (no constraints)
 * over the vertices or cells of a finite-element mesh.
 *
 * Special case of Field with a uniform fiber dimension and no constraints.
 */

#if !defined(pylith_topology_fielduniform_hh)
#define pylith_topology_fielduniform_hh

// Include directives ---------------------------------------------------
#include "Field.hh" // ISA Field

// Forward declarations -------------------------------------------------
namespace pylith {
  namespace topology {
    class FieldUniform;
    class TestFieldUniform;
  } // topology
} // pylith

// FieldUniform ---------------------------------------------------------
class pylith::topology::FieldUniform : public Field
{ // FieldUniform
  friend class TestFieldUniform; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Sieve mesh.
   */
  FieldUniform(const ALE::Obj<SieveMesh>& mesh,
	       const int fiberDim);

  /// Destructor.
  ~FieldUniform(void);

  /** Create section given points.
   *
   * @param points Mesh points over which to define section.
   */
  void createSection(const ALE::Obj<SieveMesh::label_sequence>& points);

  /** Create section given chart. This allows a chart to be reused
   * across multiple fields, reducing memory usage.
   *
   * @param chart Chart defining points over which section is defined.
   */
  void createSection(const SieveRealSection::chart_type& chart);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const int _fiberDim; ///< Fiber dimension

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  FieldUniform(const FieldUniform&); ///< Not implemented
  const FieldUniform& operator=(const FieldUniform&); ///< Not implemented

}; // FieldUniform

#endif // pylith_topology_fielduniform_hh


// End of file 
