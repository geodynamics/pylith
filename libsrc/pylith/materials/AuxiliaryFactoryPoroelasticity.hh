// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/AuxiliaryFactoryPoroelasticity.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the poroelasticity equation.
 */

#if !defined(pylith_materials_auxiliaryfactoryporoelasticity_hh)
#define pylith_materials_auxiliaryfactoryporoelasticity_hh

#include "materialsfwd.hh" // forward declarations
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
    void addGravityField(spatialdata::spatialdb::GravityField *gf);

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

    /// Add constant pressure sourcesubfield to auxiliary fields.
    void addConstantPressureSource(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPoroelasticity(const AuxiliaryFactoryPoroelasticity &); ///< Not implemented.
    const AuxiliaryFactoryPoroelasticity &operator=(const AuxiliaryFactoryPoroelasticity &); ///< Not implemented

}; // class AuxiliaryFactoryPoroelasticity

#endif // pylith_materials_auxiliaryfactoryporoelasticity_hh

// End of file
