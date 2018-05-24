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

/** @file libsrc/faults/AuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for faults.
 */

#if !defined(pylith_faults_auxiliaryfactory_hh)
#define pylith_faults_auxiliaryfactory_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

// AuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for faults.
class pylith::faults::AuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactory(void);

    /// Destructor.
    ~AuxiliaryFactory(void);

    /// Add fault strike direction subfield to auxiliary field.
    void strikeDir(void);

    /// Add fault up-dip direction modulus subfield to auxiliary field.
    void upDipDir(void);

    /// Add fault normal direction subfield to auxiliary field.
    void normalDir(void);

    /// Add slip subfield to auxiliary field.
    void slip(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    AuxiliaryFactory(const AuxiliaryFactory &);   ///< Not implemented.
    const AuxiliaryFactory& operator=(const AuxiliaryFactory&);   ///< Not implemented

}; // class AuxiliaryFactory

#endif // pylith_faults_auxiliaryfactory_hh


// End of file
