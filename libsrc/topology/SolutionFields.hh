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
 * @file libsrc/topology/SolutionFields.hh
 *
 * @brief Object for managing solution fields over a finite-element
 * mesh.
 */

#if !defined(pylith_topology_solutionfields_hh)
#define pylith_topology_solutionfields_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include <vector> // HASA std::vector

#include "Mesh.hh" // ISA Field<Mesh>
#include "Field.hh" // ISA Fields< Field<Mesh> >
#include "Fields.hh" // ISA Fields< Field<Mesh> >

// SolutionFields -------------------------------------------------------
class pylith::topology::SolutionFields : public Fields<Field<Mesh> >
{ // SolutionFields
  friend class TestSolutionFields; // unit testing

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param mesh Finite-element mesh.
   */
  SolutionFields(const Mesh& mesh);

  /// Destructor.
  ~SolutionFields(void);

  /** Set name of solution field.
   *
   * @param name Name of field that is the solution.
   */
  void solutionField(const char* name);

  /** Get solution field.
   *
   * @returns Solution field.
   */
  const Field<Mesh>& getSolution(void) const;

  /** Get solution field.
   *
   * @returns Solution field.
   */
  Field<Mesh>& getSolution(void);

  /** Create history manager for a subset of the managed fields.
   *
   * @param fields Fields in history (first is most recent).
   * @param size Number of fields in history.
   */
  void createHistory(const char** fields,
		     const int size);

  /** Shift fields in history. Handles to fields are shifted so that
   *  the most recent values become associated with the second most
   *  recent item in the history, etc.
   */
  void shiftHistory(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Name of field that corresponds to the solution.
  std::string _solutionName;

  /// History manager for a subset of the fields
  std::vector<std::string> _history;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SolutionFields(const SolutionFields&); ///< Not implemented
  const SolutionFields& operator=(const SolutionFields&); ///< Not implemented

}; // SolutionFields

#endif // pylith_topology_solutionfields_hh


// End of file 
