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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
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
    class FaultCohesive;
    class FaultCohesiveLagrange;
    class FaultCohesiveKin;
    class FaultCohesiveDyn;
    class FaultCohesiveImpulses;
    class FaultCohesiveTract;

    class EqKinSrc;
    class SlipTimeFn;
    class BruneSlipFn;
    class ConstRateSlipFn;
    class LiuCosSlipFn;
    class StepSlipFn;
    class TimeHistorySlipFn;

    class Nucleator;
    class TractPerturbation;

    class TopologyOps;
  } // faults
} // pylith


#endif // pylith_faults_bcfwd_hh


// End of file 
