// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/sources/AuxiliaryFactoryPointForce.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the wellbore source equation.
 */

#if !defined(pylith_sources_auxiliaryfactorypointforce_hh)
#define pylith_sources_auxiliaryfactorypointforce_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::sources::AuxiliaryFactoryPointForce : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryPointForce; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryPointForce(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryPointForce(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addMomentTensor(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addTimeDelay(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPointForce(const AuxiliaryFactoryPointForce &); ///< Not implemented.
    const AuxiliaryFactoryPointForce& operator=(const AuxiliaryFactoryPointForce&); ///< Not implemented

}; // class AuxiliaryFactoryPointForce

#endif // pylith_sources_auxiliaryfactorypointforce_hh

// End of file
