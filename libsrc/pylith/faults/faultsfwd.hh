// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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
        class FaultCohesive;
        class FaultCohesiveKin;
        class FaultCohesiveImpulses;
        class AuxiliaryFieldFactory;
        class DiagnosticFieldFactory;
        class DerivedFieldFactory;

        class KinSrc;
        class KinSrcConstRate;
        class KinSrcStep;
        class KinSrcRamp;
        class KinSrcBrune;
        class KinSrcLiuCos;
        class KinSrcTimeHistory;
        class KinSrcAuxiliaryFactory;

        class TopologyOps;
        class FaultOps;
    } // faults
} // pylith

#endif // pylith_faults_faultsfwd_hh

// End of file
