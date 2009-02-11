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
 * @file libsrc/topology/Visitor.hh
 *
 * @brief Visitor for managing values extracted from restriction of a
 * field to a cell.
 *
 * Wraps a Sieve RestrictVisitor.
 */

#if !defined(pylith_topology_restrictvisitor_hh)
#define pylith_topology_restrictvisitor_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// RestrictVisitor ------------------------------------------------------
template<typename field_type>
class pylith::topology::RestrictVisitor
{ // Field
  friend class TestRestrictVisitor; // unit testing

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

  // Convenience typedefs
  typedef typename field_type::RealSection RealSection;

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Constructor with field and size.
   *
   * @param field Field over finite-element mesh.
   * @param size Fiber dimension for field.
   */
  RestrictVisitor(const field_type& field,
		  const int size);

  /** Constructor with field and array.
   *
   * @param field Field over finite-element mesh.
   * @param values Array of values to use for storage.
   */
  RestrictVisitor(const field_type& field,
		  const double_array& values);

  /// Destructor.
  ~RestrictVisitor(void);

  /** Visit field for cell and retrieve values.
   *
   * @param cell Cell in finite-element mesh.
   */
  void visit(const typename field_type::Mesh::SieveMesh::point_type& cell);

  /** Get field values previously retrieved for cell.
   *
   * @returns Values for field.
   */
  const double* values(void) const;

  /// Clear values associated with cell.
  void clear(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  const field_type& _field; ///< Field associated with visitor

  /// Sieve visitor that manages storage.
  ALE::ISieveVisitor::RestrictVisitor<RealSection> _visitor;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  RestrictVisitor(const RestrictVisitor&); ///< Not implemented
  const RestrictVisitor& operator=(const RestrictVisitor&); ///< Not implemented

}; // Field

#include "RestrictVisitor.icc"

#endif // pylith_topology_restrictvisitor_hh


// End of file 
