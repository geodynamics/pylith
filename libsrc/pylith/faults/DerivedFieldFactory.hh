// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/faults/DerivedFieldFactory.hh
 *
 * @brief C++ helper class for setting up derived subfields for  faults.
 */

#if !defined(pylith_faults_derivedfieldfactory_hh)
#define pylith_faults_derivedfieldfactory_hh

#include "pylith/faults/faultsfwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA FieldFactory

class pylith::faults::DerivedFieldFactory : public pylith::topology::FieldFactory {
    friend class TestDerivedFieldFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFieldFactory(void);

    /// Destructor.
    virtual ~DerivedFieldFactory(void);

    /// Add traction change subfield to derived field.
    void addTractionChange(void);

    /// Add subfields using discretizations provided.
    virtual void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DerivedFieldFactory(const DerivedFieldFactory &); ///< Not implemented.
    const DerivedFieldFactory& operator=(const DerivedFieldFactory&); ///< Not implemented

}; // class DerivedFieldFactory

#endif // pylith_faults_derivedfieldfactory_hh

// End of file
