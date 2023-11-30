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

/** @file libsrc/bc/DiagnosticFieldFactory.hh
 *
 * @brief C++ helper class for setting up derived subfields for bc.
 */

#if !defined(pylith_bc_diagnosticfieldfactory_hh)
#define pylith_bc_diagnosticfieldfactory_hh

#include "bcfwd.hh" // forward declarations
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

#endif // pylith_bc_diagnosticfieldfactory_hh

// End of file
