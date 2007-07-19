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

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

#include "Fault.hh" // ISA Fault

#include "pylith/utils/sievefwd.hh" // HOLDSA PETSc Mesh
#include "pylith/utils/petscfwd.h" // USES PetscMat

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesive;
    class TestFaultCohesive; // unit testing
  } // faults
} // pylith

/// C++ abstract base class for a fault surface implemented with
/// cohesive elements.
class pylith::faults::FaultCohesive : public Fault
{ // class FaultCohesive
  friend class TestFaultCohesive; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesive(void);

  /// Destructor.
  virtual
  ~FaultCohesive(void);

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   */
  void adjustTopology(const ALE::Obj<ALE::Mesh>& mesh);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  virtual
  bool _useLagrangeConstraints(void) const = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesive(const FaultCohesive& m);

  /// Not implemented
  const FaultCohesive& operator=(const FaultCohesive& m);

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  ALE::Obj<ALE::Mesh>* _faultMesh;

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh


// End of file 
