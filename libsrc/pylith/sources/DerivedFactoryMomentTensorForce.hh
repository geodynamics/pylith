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

/** @file libsrc/sources/DerivedFactoryMomentTensorForce.hh
 *
 * @brief C++ helper class for setting up derived subfields for elastic sources.
 */

#if !defined(pylith_sources_derivedfactorymomenttensorforce_hh)
#define pylith_sources_derivedfactorymomenttensorforce_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA AuxiliaryFactory

class pylith::sources::DerivedFactoryMomentTensorForce : public pylith::topology::FieldFactory {
    friend class TestDerivedFactoryMomentTensorForce; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFactoryMomentTensorForce(void);

    /// Destructor.
    virtual ~DerivedFactoryMomentTensorForce(void);

    /// Add Cauchy stress subfield to derived field.
    void addCauchyStress(void);

    /// Add Cauchy (infinitesimal) strain subfield to derived field.
    void addCauchyStrain(void);

    /// Add subfields using discretizations provided.
    void addSubfields(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    DerivedFactoryMomentTensorForce(const DerivedFactoryMomentTensorForce &); ///< Not implemented.
    const DerivedFactoryMomentTensorForce& operator=(const DerivedFactoryMomentTensorForce&); ///< Not implemented

}; // class DerivedFactoryMomentTensorForce

#endif // pylith_materials_derivedfactorymomenttensorforce_hh

// End of file
