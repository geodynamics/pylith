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
 * @file pylith/meshio/CellFilter.hh
 *
 * @brief C++ object for filtering cell fields when outputing
 * finite-element data.
 */

#if !defined(pylith_meshio_cellfilter_hh)
#define pylith_meshio_cellfilter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/feassemble/topologyfwd.hh" // HOLDSA Quadrature<Mesh>
#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Quadrature<Mesh>

// CellFilter -----------------------------------------------------------
template<typename mesh_type>
class pylith::meshio::CellFilter
{ // CellFilter

// PROTECTED TYPEDEFS ///////////////////////////////////////////////////
protected:

  // Convenience typedefs
  typedef typename mesh_type::RealSection RealSection;
  typedef typename mesh_type::SieveMesh SieveMesh;
  typedef typename SieveMesh::label_sequence label_sequence;
  typedef typename RealSection::chart_type chart_type;

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  CellFilter(void);

  /// Destructor
  ~CellFilter(void);

  /** Create copy of filter.
   *
   * @returns Copy of filter.
   */
  virtual
  CellFilter* clone(void) const = 0;

  /** Set quadrature associated with cells.
   *
   * @param q Quadrature for cells.
   */
  void quadrature(const feassemble::Quadrature<mesh_type>* q);

  /** Filter field. Field type of filtered field is returned via an argument.
   *
   * @param fieldIn Field to filter.
   * @param label Value of label of cells to filter.
   * @param labelId Id associated with label of cells to filter.
   *
   * @returns Averaged field.
   */
  virtual
  const topology::Field<mesh_type>&
  filter(const topology::Field<mesh_type>& fieldIn,
	 const char* label =0,
	 const int labelId =0) = 0;

// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f Filter to copy.
   * @returns Pointer to this.
   */
  CellFilter(const CellFilter& f);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  const CellFilter& operator=(const CellFilter&);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Quadrature associated with cells.
  feassemble::Quadrature<mesh_type>* _quadrature;

}; // CellFilter

#include "CellFilter.cc" // template definitions

#endif // pylith_meshio_cellfilter_hh


// End of file 
