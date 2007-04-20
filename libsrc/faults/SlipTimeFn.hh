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

/** @file libsrc/faults/SlipTimeFn.hh
 *
 * @brief C++ abstract base class for kinematic slip time function.
 *
 * Interface definition for slip time function.
 */

#if !defined(pylith_faults_sliptimefn_hh)
#define pylith_faults_sliptimefn_hh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class SlipTimeFn;
    class TestSlipTimeFn; // unit testing
  } // faults
} // pylith

/// C++ abstract base class for Fault object.
class pylith::faults::SlipTimeFn
{ // class SlipTimeFn
  friend class TestSlipTimeFn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  SlipTimeFn(void);

  /// Destructor.
  virtual
  ~SlipTimeFn(void);

  /** Create copy of fault.
   *
   * @returns Copy of fault.
   */
  virtual
  SlipTimeFn* clone(void) const = 0;

  /** Compute slip using slip time function.
   *
   * @param t Time relative to slip starting time at point
   * @param finalSlip Final slip at point
   * @param peakRate Peak slip rate at point
   *
   * @returns Slip at point at time t
   */
  virtual
  double compute(const double t,
		 const double finalSlip,
		 const double peakRate) const = 0;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param f SlipTimeFn to copy
   */
  SlipTimeFn(const SlipTimeFn& f);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const SlipTimeFn& operator=(const SlipTimeFn& f);

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
