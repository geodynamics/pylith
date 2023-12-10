// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::materials::AuxiliaryFactoryPoroelasticity : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryPoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryPoroelasticity(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryPoroelasticity(void);

    /// Add body force subfield to auxiliary subfields.
    void addBodyForce(void);

    /** Add gravity subfield to auxiliary subfields.
     *
     * @param[in] gf Gravity field.
     */
    void addGravityField(spatialdata::spatialdb::GravityField* gf);

    /// Add porosity subfield to auxiliary subfields.
    void addPorosity(void);

    /// Add solid density subfield to auxiliary subfields.
    void addSolidDensity(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensity(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addFluidViscosity(void);

    /// Add reference sourceDensity subfield to auxiliary fields.
    void addSourceDensity(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPoroelasticity(const AuxiliaryFactoryPoroelasticity &); ///< Not implemented.
    const AuxiliaryFactoryPoroelasticity& operator=(const AuxiliaryFactoryPoroelasticity&); ///< Not implemented

}; // class AuxiliaryFactoryPoroelasticity

// End of file
