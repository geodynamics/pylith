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

/** @file libsrc/faults/faultsfwd.hh
 *
 * @brief Forward declarations for PyLith faults objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_faults_faultsfwd_hh)
#define pylith_faults_faultsfwd_hh

namespace pylith {
  namespace faults {

    class CohesiveTopology;

    class Fault;
    class FaultCohesiveDyn;
    class FaultCohesiveKin;

    class EqKinSrc;
    class SlipTimeFn;
    class BruneSlipFn;
    class ConstRateSlipFn;
    class LiuCosSlipFn;
    class StepSlipFn;

    class TopologyOps;
    template<typename Sieve, typename Renumbering> class ReplaceVisitor;
    template<typename Sieve> class ClassifyVisitor;

  } // faults
} // pylith


#endif // pylith_faults_bcfwd_hh


// End of file 
