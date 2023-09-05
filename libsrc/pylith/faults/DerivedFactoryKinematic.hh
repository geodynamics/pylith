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

/** @file libsrc/faults/DerivedFactoryKinematic.hh
 *
 * @brief C++ helper class for setting up derived subfields for kinematic faults.
 */

#if !defined(pylith_faults_derivedfactorykinematic_hh)
#define pylith_faults_derivedfactorykinematic_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::DerivedFactoryKinematic : public pylith::topology::FieldFactory {
    friend class TestDerivedFactoryKinematic; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFactoryKinematic(void);

    /// Destructor.
    virtual ~DerivedFactoryKinematic(void);

    /// Add traction change subfield to derived field.
    void addTractionChange(void);

    /// Add subfields using discretizations provided.
    virtual void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DerivedFactoryKinematic(const DerivedFactoryKinematic &); ///< Not implemented.
    const DerivedFactoryKinematic& operator=(const DerivedFactoryKinematic&); ///< Not implemented

}; // class DerivedFactoryKinematic

#endif // pylith_faults_derivedfactorykinematic_hh

// End of file
