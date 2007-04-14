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

/** @file libsrc/faults/FaultCohesiveKin.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

#include "FaultCohesive.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class FaultCohesiveKin;
    class TestFaultCohesiveKin; // unit testing
  } // faults
} // pylith

/// C++ implementation for a fault surface with kinematic (prescribed)
/// slip implemented with cohesive elements.
class pylith::faults::FaultCohesiveKin : public FaultCohesive
{ // class FaultCohesiveKin
  friend class TestFaultCohesiveKin; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveKin(void);

  /// Destructor.
  virtual
  ~FaultCohesiveKin(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  Fault* clone(void) const;  

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Fault to copy
   */
  FaultCohesiveKin(const FaultCohesiveKin& m);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

}; // class FaultCohesiveKin

#include "FaultCohesiveKin.icc" // inline methods

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
