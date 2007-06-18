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

/** @file libsrc/faults/Fault.hh
 *
 * @brief C++ abstract base class for Fault object.
 *
 * Interface definition for fault.
 */

#if !defined(pylith_faults_fault_hh)
#define pylith_faults_fault_hh

#include "pylith/utils/sievefwd.hh" // USES PETSc Mesh
#include "pylith/utils/arrayfwd.hh" // USES double_array

#include <string> // HASA std::string

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class Fault;
    class TestFault; // unit testing
  } // faults
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

/// C++ abstract base class for Fault object.
class pylith::faults::Fault
{ // class Fault
  friend class TestFault; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Fault(void);

  /// Destructor.
  virtual
  ~Fault(void);

  /** Set identifier of fault.
   *
   * @param value Fault identifier
   */
  void id(const int value);

  /** Get identifier of fault.
   *
   * @returns Fault identifier
   */
  int id(void) const;

  /** Set label of fault.
   *
   * @param value Label of fault
   */
  void label(const char* value);

  /** Get label of fault.
   *
   * @returns Label of fault
   */
  const std::string& label(void) const;

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   */
  virtual
  void adjustTopology(const ALE::Obj<ALE::Mesh>& mesh) = 0;

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces in a 3-D domain).
   * @param normalDir General preferred direction for fault normal
   *   (used to pick which of two possible normal directions for
   *   interface; only applies to fault surfaces in a 3-D domain).
   */
  virtual
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir,
		  const double_array& normalDir) = 0;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Fault(const Fault& m);

  /// Not implemented
  const Fault& operator=(const Fault& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  int _id; ///< Fault identifier
  std::string _label; ///< Label of fault

}; // class Fault

#include "Fault.icc" // inline methods

#endif // pylith_faults_fault_hh


// End of file 
