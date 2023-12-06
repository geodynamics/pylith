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

/** @file libsrc/faults/AuxiliaryFieldFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for faults.
 */

#if !defined(pylith_faults_auxiliaryfactorykinematic_hh)
#define pylith_faults_auxiliaryfactorykinematic_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::AuxiliaryFieldFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFieldFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFieldFactory(void);

    /// Destructor.
    ~AuxiliaryFieldFactory(void);

    /// Add slip subfield to auxiliary field.
    void addSlip(void);

    /// Add slip rate subfield to auxiliary field.
    void addSlipRate(void);

    /// Add slip acceleration subfield to auxiliary field.
    void addSlipAcceleration(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFieldFactory(const AuxiliaryFieldFactory &); ///< Not implemented.
    const AuxiliaryFieldFactory& operator=(const AuxiliaryFieldFactory&); ///< Not implemented

}; // class AuxiliaryFieldFactory

#endif // pylith_faults_auxiliaryfactorykinematic_hh

// End of file
