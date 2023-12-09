// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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
