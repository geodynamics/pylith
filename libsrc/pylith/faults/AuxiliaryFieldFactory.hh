// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
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
