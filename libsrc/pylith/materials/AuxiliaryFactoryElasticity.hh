// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::materials::AuxiliaryFactoryElasticity : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryElasticity(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryElasticity(void);

    /// Add density subfield to auxiliary subfields.
    void addDensity(void);

    /// Add body force subfield to auxiliary subfields.
    void addBodyForce(void);

    /** Add gravity subfield to auxiliary subfields.
     *
     * @param[in] gf Gravity field.
     */
    void addGravityField(spatialdata::spatialdb::GravityField* gf);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryElasticity(const AuxiliaryFactoryElasticity &); ///< Not implemented.
    const AuxiliaryFactoryElasticity& operator=(const AuxiliaryFactoryElasticity&); ///< Not implemented

}; // class AuxiliaryFactoryElasticity

// End of file
