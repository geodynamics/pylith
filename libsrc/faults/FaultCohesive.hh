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

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
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
   * If firstFaultVertex == 0, then firstFaultVertex is set to the first point
   * not currently used in the mesh, and firstFaultCell is incremented with this
   * point. These values are updated as new fault vertices and cells are added.
   *
   * @param mesh PETSc mesh.
   * @param firstFaultVertex The first point eligible to become a new fault vertex
   * @param firstFaultCell The first point eligible to become a new fault cell
   * @param flipFault Flip fault orientation.
   */
  void adjustTopology(topology::Mesh* const mesh,
                      int *firstFaultVertex,
                      int *firstFaultCell,
                      const bool flipFault = false);

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  virtual
  bool useLagrangeConstraints(void) const = 0;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// If true, use fault mesh to define fault; otherwise, use group of
  /// vertices to define fault.
  bool _useFaultMesh;

  std::string _faultMeshFilename; /// Filename for fault mesh UCD file.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  FaultCohesive(const FaultCohesive&); ///< Not implemented
  const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh


// End of file 
