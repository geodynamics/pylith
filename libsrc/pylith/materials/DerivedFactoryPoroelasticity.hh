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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/DerivedFactoryElasticity.hh
 *
 * @brief C++ helper class for setting up derived subfields for elastic materials.
 */

#if !defined(pylith_materials_derivedfactoryporoelasticity_hh)
#define pylith_materials_derivedfactoryporoelasticity_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/materials/DerivedFactoryElasticity.hh" // ISA DerivedFactoryElasticity

class pylith::materials::DerivedFactoryPoroelasticity : public pylith::materials::DerivedFactoryElasticity {
    friend class TestDerivedFactoryElasticity; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFactoryPoroelasticity(void);

    /// Destructor.
    virtual ~DerivedFactoryPoroelasticity(void);

    /// Add bulk density subfield to derived field
    void addBulkDensity(void);

    /// Add water content subfield to derived field
    void addWaterContent(void);

    /// Add subfields using discretizations provided.
    void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DerivedFactoryPoroelasticity(const DerivedFactoryPoroelasticity &); ///< Not implemented.
    const DerivedFactoryPoroelasticity& operator=(const DerivedFactoryPoroelasticity&); ///< Not implemented

}; // class DerivedFactoryPoroelasticity

#endif // pylith_materials_derivedfactoryporoelasticity_hh

// End of file
