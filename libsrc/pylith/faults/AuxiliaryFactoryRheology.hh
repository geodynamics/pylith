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

/** @file libsrc/faults/AuxiliaryFactoryRheology.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the fault constitutive models.
 */

#if !defined(pylith_faults_auxiliaryfactoryelasticity_hh)
#define pylith_faults_auxiliaryfactoryelasticity_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::AuxiliaryFactoryRheology : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryRheology; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryRheology(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryRheology(void);

    /// Add cohesion to auxiliary subfields.
    void addCohesion(void);

    /// Add static friction coefficient to auxiliary subfields.
    void addStaticCoefficient(void);

    /// Add dynamic friction coefficient to auxiliary subfields.
    void addDynamicCoefficient(void);

    /// Add slip-weakening parameter to auxiliary subfields.
    void addSlipWeakeningParameter(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryRheology(const AuxiliaryFactoryRheology &); ///< Not implemented.
    const AuxiliaryFactoryRheology& operator=(const AuxiliaryFactoryRheology&); ///< Not implemented

}; // class AuxiliaryFactoryRheology

#endif // pylith_faults_auxiliaryfactoryelasticity_hh

// End of file
