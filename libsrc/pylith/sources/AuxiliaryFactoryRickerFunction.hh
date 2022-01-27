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

/** @file libsrc/sources/AuxiliaryFactoryRickerFunction.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the wellbore source equation.
 */

#if !defined(pylith_sources_auxiliaryfactoryrickerfunction_hh)
#define pylith_sources_auxiliaryfactoryrickerfunction_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/sources/AuxiliaryFactoryPointForce.hh" // ISA AuxiliaryFactory

class pylith::sources::AuxiliaryFactoryRickerFunction : public pylith::sources::AuxiliaryFactoryPointForce {
    friend class TestAuxiliaryFactoryRickerFunction; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryRickerFunction(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryRickerFunction(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addRickerCenterFrequency(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryRickerFunction(const AuxiliaryFactoryRickerFunction &); ///< Not implemented.
    const AuxiliaryFactoryRickerFunction& operator=(const AuxiliaryFactoryRickerFunction&); ///< Not implemented

}; // class AuxiliaryFactoryRickerFunction

#endif // pylith_sources_auxiliaryfactoryrickerfunction_hh

// End of file
