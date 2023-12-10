// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/bc/bcfwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA AuxiliaryFactory

class pylith::bc::DiagnosticFieldFactory : public pylith::topology::FieldFactory {
    friend class TestDerivedFieldFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DiagnosticFieldFactory(void);

    /// Destructor.
    virtual ~DiagnosticFieldFactory(void);

    /// Add boundary normal direction subfield to diagnostic field.
    void addNormalDir(void);

    /// Add (horizontla) tangential direction subfield to diagnostic field.
    void addTangentialDirHoriz(void);

    /// Add (vertical) tangential direction subfield to diagnostic field.
    void addTangentialDirVert(void);

    /// Add subfields using discretizations provided.
    virtual void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DiagnosticFieldFactory(const DiagnosticFieldFactory &); ///< Not implemented.
    const DiagnosticFieldFactory& operator=(const DiagnosticFieldFactory&); ///< Not implemented

}; // class DiagnosticFieldFactory

// End of file
