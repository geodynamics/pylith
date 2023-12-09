// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/faults/DiagnosticFieldFactory.hh
 *
 * @brief C++ helper class for setting up derived subfields for faults.
 */

#if !defined(pylith_faults_diagnosticfieldfactory_hh)
#define pylith_faults_diagnosticfieldfactory_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::DiagnosticFieldFactory : public pylith::topology::FieldFactory {
    friend class TestDerivedFieldFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DiagnosticFieldFactory(void);

    /// Destructor.
    virtual ~DiagnosticFieldFactory(void);

    /// Add fault normal direction subfield to auxiliary field.
    void addNormalDir(void);

    /// Add fault strike direction subfield to auxiliary field.
    void addStrikeDir(void);

    /// Add fault up-dip direction modulus subfield to auxiliary field.
    void addUpDipDir(void);

    /// Add subfields using discretizations provided.
    virtual void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DiagnosticFieldFactory(const DiagnosticFieldFactory &); ///< Not implemented.
    const DiagnosticFieldFactory& operator=(const DiagnosticFieldFactory&); ///< Not implemented

}; // class DiagnosticFieldFactory

#endif // pylith_faults_diagnosticfieldfactory_hh

// End of file
