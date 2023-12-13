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

/** @file libsrc/sources/AuxiliaryFactoryWellboreSource.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the wellbore source equation.
 */

#if !defined(pylith_sources_auxiliaryfactorywellboresource_hh)
#define pylith_sources_auxiliaryfactorywellboresource_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::sources::AuxiliaryFactoryWellboreSource : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryWellboreSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryWellboreSource(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryWellboreSource(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensity(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addFluidViscosity(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add wellbore radius subfield to auxiliary subfields.
    void addWellboreRadius(void);

    /// Add wellbore length subfield to auxiliary subfields.
    void addWellboreLength(void);

    /// Add wellbore pressure subfield to auxiliary subfields.
    void addWellborePressure(void);

    /// Add wellbore character subfield to auxiliary subfields.
    void addWellboreCharacter(void);

    /// Add element dimensions subfield to auxiliary subfields.
    void addElementDimensions(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryWellboreSource(const AuxiliaryFactoryWellboreSource &); ///< Not implemented.
    const AuxiliaryFactoryWellboreSource& operator=(const AuxiliaryFactoryWellboreSource&); ///< Not implemented

}; // class AuxiliaryFactoryWellboreSource

#endif // pylith_sources_auxiliaryfactorywellboresource_hh

// End of file
