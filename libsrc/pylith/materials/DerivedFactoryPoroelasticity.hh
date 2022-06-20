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

/** @file libsrc/materials/DerivedFactoryPoroelasticity.hh
 *
 * @brief C++ helper class for setting up derived subfields for elastic materials.
 */

#if !defined(pylith_materials_derivedfactoryporoelasticity_hh)
#define pylith_materials_derivedfactoryporoelasticity_hh

#include "materialsfwd.hh"                 // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA AuxiliaryFactory

class pylith::materials::DerivedFactoryPoroelasticity : public pylith::topology::FieldFactory
{
    friend class TestDerivedFactoryPoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    DerivedFactoryPoroelasticity(void);

    /// Destructor.
    virtual ~DerivedFactoryPoroelasticity(void);

    /// Add Cauchy stress subfield to derived field.
    void addCauchyStress(void);

    /// Add Cauchy (infinitesimal) strain subfield to derived field.
    void addCauchyStrain(void);

    /// Add add generated porosity subfield to derived field.
    void addOutputPorosity(void);

    /// Add subfields using discretizations provided.
    void addSubfields(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    DerivedFactoryPoroelasticity(const DerivedFactoryPoroelasticity &);                  ///< Not implemented.
    const DerivedFactoryPoroelasticity &operator=(const DerivedFactoryPoroelasticity &); ///< Not implemented

}; // class DerivedFactoryPoroelasticity

#endif // pylith_materials_derivedfactoryporoelasticity_hh

// End of file
