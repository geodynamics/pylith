// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, Rice University
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
 * @file libsrc/feassemble/IntegratorPointwise.hh
 *
 * @brief Object containing operations for implicit and explicit
 * time integration of the equations defined by pointwise functions.
 */

#if !defined(pylith_feassemble_integratorpointwise_hh)
#define pylith_feassemble_integratorpointwise_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "Integrator.hh" // ISA Integrator

#include "pylith/topology/topologyfwd.hh" // HOLDSA Field

// IntegratorPointwise -------------------------------------------------
/** @brief General operations for implicit and explicit
 * time integration of equations defined by pointwise functions.
 */
class pylith::feassemble::IntegratorPointwise : public Integrator
{ // IntegratorPointwise
  friend class TestIntegratorPointwise; // unit testing

// PUBLIC TYPEDEFS //////////////////////////////////////////////////////
public :

// PUBLIC MEMBERS ///////////////////////////////////////////////////////
public :

  /// Constructor
  IntegratorPointwise(void);

  /// Destructor
  virtual
  ~IntegratorPointwise(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);

  /** Get auxiliary fields.
   *
   * @return field Field over material.
   */
  const topology::Field& auxFields(void) const;

  /** Check whether material has a given auxiliary field.
   *
   * @param name Name of field.
   *
   * @returns True if material has auxiliary field, false otherwise.
   */
  bool hasAuxField(const char* name);

  /** Get auxiliary field.
   *
   * @param field Field over material.
   * @param name Name of field to retrieve.
   */
  void getAuxField(topology::Field *field,
		   const char* name) const;

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void updateStateVars(const PylithScalar t,
		       topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const;
  
// PROTECTED METHODS ////////////////////////////////////////////////////
protected :

  /** Set residual and Jacobian kernels.
   *
   * @param solution Solution field.
   */
  virtual
  void _setFEKernels(const topology::Field& solution) const = 0;


// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Auxiliary fields for this problem
  topology::Field *_auxFields;

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  /// Not implemented.
  IntegratorPointwise(const IntegratorPointwise&);

  /// Not implemented
  const IntegratorPointwise& operator=(const IntegratorPointwise&);

}; // IntegratorPointwise

#endif // pylith_feassemble_integratorpointwise_hh


// End of file 
