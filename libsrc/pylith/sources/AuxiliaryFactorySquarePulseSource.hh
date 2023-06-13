// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/sources/AuxiliaryFactorySquarePulseSource.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the squarepulse source equation.
 */

#if !defined(pylith_sources_auxiliaryfactorysquarepulsesource_hh)
#define pylith_sources_auxiliaryfactorysquarepulsesource_hh

#include "sourcesfwd.hh"                         // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::sources::AuxiliaryFactorySquarePulseSource : public pylith::feassemble::AuxiliaryFactory
{
    friend class TestAuxiliaryFactorySquarePulseSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    AuxiliaryFactorySquarePulseSource(void);

    /// Destructor.
    virtual ~AuxiliaryFactorySquarePulseSource(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addVolumeFlowRate(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    AuxiliaryFactorySquarePulseSource(const AuxiliaryFactorySquarePulseSource &);                  ///< Not implemented.
    const AuxiliaryFactorySquarePulseSource &operator=(const AuxiliaryFactorySquarePulseSource &); ///< Not implemented

}; // class AuxiliaryFactorySquarePulseSource

#endif // pylith_sources_auxiliaryfactorysquarepulsesource_hh

// End of file
