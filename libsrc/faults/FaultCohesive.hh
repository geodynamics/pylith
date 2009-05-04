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

// Include directives ---------------------------------------------------
#include "Fault.hh" // ISA Fault

// FaultCohesive --------------------------------------------------------
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

  /** Set flag for using fault mesh or group of vertices to define
   * fault surface.
   *
   * @param flag True if using fault mesh, false if using vertices.
   */
  void useFaultMesh(const bool flag);

  // TEMPORARY
  /** Set filename of UCD file for fault mesh.
   *
   * @param filename Filename for UCD file.
   */
  void faultMeshFilename(const char* filename);

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh.
   * @param flipFault Flip fault orientation.
   */
  void adjustTopology(topology::Mesh* const mesh,
		      const bool flipFault =false);

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

  FaultCohesive(const FaultCohesive&); ///< Not implemented
  const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// If true, use fault mesh to define fault; otherwise, use group of
  /// vertices to define fault.
  bool _useFaultMesh;

  std::string _faultMeshFilename; /// Filename for fault mesh UCD file.

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh


// End of file 
