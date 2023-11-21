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
