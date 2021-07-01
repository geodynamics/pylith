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

/** @file libsrc/materials/AuxiliaryFactoryElastic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for elastic materials.
 */

#if !defined(pylith_materials_auxiliaryfactoryelastic_hh)
#define pylith_materials_auxiliaryfactoryelastic_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // ISA AuxiliaryFactoryElasticity

class pylith::materials::AuxiliaryFactoryElastic : public pylith::materials::AuxiliaryFactoryElasticity {
    friend class TestAuxiliaryFactoryElastic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryElastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryElastic(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add bulk subfield to auxiliary subfields.
    void addBulkModulus(void);

    /// Add reference stress subfield to auxiliary subfields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary subfields.
    void addReferenceStrain(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryElastic(const AuxiliaryFactoryElastic &); ///< Not implemented.
    const AuxiliaryFactoryElastic& operator=(const AuxiliaryFactoryElastic&); ///< Not implemented

}; // class AuxiliaryFactoryElastic

#endif // pylith_materials_auxiliaryfactoryelastic_hh

// End of file
