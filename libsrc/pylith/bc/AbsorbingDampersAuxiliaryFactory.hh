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

/** @file libsrc/bc/AbsorbingDampersAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for absorbing boundary conditions.
 */

#if !defined(pylith_bc_absorbingdampersauxiliaryfactory_hh)
#define pylith_bc_absorbingdampersauxiliaryfactory_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

// AbsorbingDampersAuxiliaryFactory-----------------------------------------------
/// @brief C++ helper class for setting up auxiliary fields for time-dependent boundary conditions.
class pylith::bc::AbsorbingDampersAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestDirichletAuxiliaryFactory;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    AbsorbingDampersAuxiliaryFactory(void);

    /// Destructor.
    ~AbsorbingDampersAuxiliaryFactory(void);

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Add density field to auxiliary fields.
    void density(void);

    /// Add shear wave speed field to auxiliary fields.
    void vs(void);

    /// Add dilatational wave speed field to auxiliary fields.
    void vp(void);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    static const char* _genericComponent; ///< Name of generic component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    AbsorbingDampersAuxiliaryFactory(const AbsorbingDampersAuxiliaryFactory &);   ///< Not implemented.
    const AbsorbingDampersAuxiliaryFactory& operator=(const AbsorbingDampersAuxiliaryFactory&);   ///< Not implemented

}; // class AbsorbingDampersAuxiliaryFactory

#endif // pylith_bc_absorbingdampersauxiliaryfactory_hh


// End of file
