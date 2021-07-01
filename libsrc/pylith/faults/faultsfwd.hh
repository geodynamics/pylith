// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
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
        class FaultCohesive;
        class FaultCohesiveKin;
        class AuxiliaryFactoryKinematic;

        class KinSrc;
        class KinSrcConstRate;
        class KinSrcStep;
        class KinSrcRamp;
        class KinSrcBrune;
        class KinSrcLiuCos;
        class KinSrcTimeHistory;
        class KinSrcAuxiliaryFactory;

        class TopologyOps;
    } // faults
} // pylith

#endif // pylith_faults_bcfwd_hh

// End of file
