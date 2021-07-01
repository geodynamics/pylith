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

/** @file libsrc/materials/AuxiliaryFactoryElasticity.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the elasticity equation.
 */

#if !defined(pylith_materials_auxiliaryfactoryelasticity_hh)
#define pylith_materials_auxiliaryfactoryelasticity_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::materials::AuxiliaryFactoryElasticity : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryElasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
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

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryElasticity(const AuxiliaryFactoryElasticity &); ///< Not implemented.
    const AuxiliaryFactoryElasticity& operator=(const AuxiliaryFactoryElasticity&); ///< Not implemented

}; // class AuxiliaryFactoryElasticity

#endif // pylith_materials_auxiliaryfactoryelasticity_hh

// End of file
