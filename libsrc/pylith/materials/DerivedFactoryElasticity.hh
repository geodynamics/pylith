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

/** @file libsrc/materials/DerivedFactoryElasticity.hh
 *
 * @brief C++ helper class for setting up derived subfields for elastic materials.
 */

#if !defined(pylith_materials_derivedfactoryelasticity_hh)
#define pylith_materials_derivedfactoryelasticity_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA AuxiliaryFactory

class pylith::materials::DerivedFactoryElasticity : public pylith::topology::FieldFactory {
    friend class TestDerivedFactoryElasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFactoryElasticity(void);

    /// Destructor.
    virtual ~DerivedFactoryElasticity(void);

    /// Add Cauchy stress subfield to derived field.
    void addCauchyStress(void);

    /// Add Cauchy (infinitesimal) strain subfield to derived field.
    void addCauchyStrain(void);

    /// Add subfields using discretizations provided.
    void addSubfields(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    DerivedFactoryElasticity(const DerivedFactoryElasticity &); ///< Not implemented.
    const DerivedFactoryElasticity& operator=(const DerivedFactoryElasticity&); ///< Not implemented

}; // class DerivedFactoryElasticity

#endif // pylith_materials_derivedfactoryelasticity_hh

// End of file
