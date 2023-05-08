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
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/AuxiliaryFactoryDynamic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for dynamic faults.
 */

#if !defined(pylith_faults_auxiliaryfactorydynamic_hh)
#define pylith_faults_auxiliaryfactorydynamic_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::AuxiliaryFactoryDynamic : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryDynamic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryDynamic(void);

    /// Destructor.
    ~AuxiliaryFactoryDynamic(void);

    /// Add traction perturbation subfield to auxiliary field.
    void addTractionPerturbation(void);

    /// Add rheology traction subfield to auxiliary field.
    void addTractionRheology(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryDynamic(const AuxiliaryFactoryDynamic &); ///< Not implemented.
    const AuxiliaryFactoryDynamic& operator=(const AuxiliaryFactoryDynamic&); ///< Not implemented

}; // class AuxiliaryFactoryDynamic

#endif // pylith_faults_auxiliaryfactorydynamic_hh

// End of file
