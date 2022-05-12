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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file libsrc/sources/sourcesfwd.hh
 *
 * @brief Forward declarations for PyLith sources objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_sources_sourcesfwd_hh)
#define pylith_sources_sourcesfwd_hh

namespace pylith {
    namespace sources {
        // New stuff
        class Source;

        class PointForce;
        class AuxiliaryFactoryPointForce;
        class DerivedFactoryPointForce;
        
        class SourceTimeFunctionPointForce;
        class AuxiliaryFactorySourceTime;
        class RickerWavelet;

        class WellboreSource;        
        class AuxiliaryFactoryWellboreSource;
        

    } // sources
} // pylith

#endif // pylith_sources_sourcesfwd_hh

// End of file
