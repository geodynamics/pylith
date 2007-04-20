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

/** @file libsrc/faults/BruneSlipFn.hh
 *
 * @brief C++ implementation of Brune slip time function.
 *
 * Slip time function follows the integral of Brune's (1907) far-field
 * time function.
 *
 * Normalize slip = 1 - exp(-t/tau)(1 + t/tau),
 * where tau = finalSlip / (exp(1.0) * peakRate)
 */

#if !defined(pylith_faults_bruneslipfn_hh)
#define pylith_faults_bruneslipfn_hh

#include "SlipTimeFn.hh"

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class BruneSlipFn;
    class TestBruneSlipFn; // unit testing
  } // faults
} // pylith

/// C++ implementation of Brune slip time function.
class pylith::faults::BruneSlipFn : public SlipTimeFn
{ // class BruneSlipFn
  friend class TestBruneSlipFn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  BruneSlipFn(void);

  /// Destructor.
  virtual
  ~BruneSlipFn(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  SlipTimeFn* clone(void) const;  

  /** Compute slip using slip time function.
   *
   * @param t Time relative to slip starting time at point
   * @param finalSlip Final slip at point
   * @param peakRate Peak slip rate at point
   *
   * @returns Slip at point at time t
   */
  double compute(const double t,
		 const double finalSlip,
		 const double peakRate) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f BruneSlipFn to copy
   */
  BruneSlipFn(const BruneSlipFn& m);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const BruneSlipFn& operator=(const BruneSlipFn& f);

}; // class FaultCohesiveKin

#include "BruneSlipFn.icc" // inline methods

#endif // pylith_faults_bruneslipfn_hh


// End of file 
