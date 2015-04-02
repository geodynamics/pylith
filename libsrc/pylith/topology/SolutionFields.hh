// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
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

#include "Fields.hh" // ISA Fields

#include <vector> // HASA std::vector

// SolutionFields -------------------------------------------------------
/** @brief Object for managing solution fields over a finite-element
 * mesh.
 */
class pylith::topology::SolutionFields : public Fields
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

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set name of field used in the solve.
   *
   * @param name Name of field that holds the solution.
   */
  void solutionName(const char* name);

  /** Get field used in solve.
   *
   * @returns Solution field.
   */
  const Field& solution(void) const;

  /** Get field used in solve.
   *
   * @returns Solution field.
   */
  Field& solution(void);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Name of field that corresponds to the "working" solution to the
  /// problem.
  std::string _solutionName;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  SolutionFields(const SolutionFields&); ///< Not implemented
  const SolutionFields& operator=(const SolutionFields&); ///< Not implemented

}; // SolutionFields

#endif // pylith_topology_solutionfields_hh


// End of file 
