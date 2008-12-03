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
 * @file pylith/topology/FieldOps.hh
 *
 * @brief Temporary object for doing operations on a PETSc
 * field. Object will be replaced by a PyLith Field object that inherits
 * or templates over the PETSc Field object.
 */

#if !defined(pylith_topology_fieldops_hh)
#define pylith_topology_fieldops_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc real_section_type

namespace pylith {
  namespace topology {
    class FieldOps;
    class TestFieldOps;
  } // topology
} // pylith

class pylith::topology::FieldOps
{ // MeshOps
  friend class TestFieldOps; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Copy values from one section to another. Sections must be
   * compatible in size and shape.
   *
   * @param dest Section to copy values into.
   * @param src Section to copy values from.
   */
  static
  void copyValues(const ALE::Obj<real_section_type>& dest,
		  const ALE::Obj<real_section_type>& src);


// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented
  FieldOps(const FieldOps&);

  /// Not implemented
  const FieldOps& operator=(const FieldOps&);


}; // FieldOps

#endif // pylith_topology_fieldops_hh


// End of file 
