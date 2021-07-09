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
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/AbsorbingDampersAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for absorbing boundary conditions.
 */

#if !defined(pylith_bc_absorbingdampersauxiliaryfactory_hh)
#define pylith_bc_absorbingdampersauxiliaryfactory_hh

#include "bcfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::bc::AbsorbingDampersAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestDirichletAuxiliaryFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AbsorbingDampersAuxiliaryFactory(void);

    /// Destructor.
    virtual ~AbsorbingDampersAuxiliaryFactory(void);

    /// Add density field to auxiliary fields.
    void addDensity(void);

    /// Add shear wave speed field to auxiliary fields.
    void addVs(void);

    /// Add dilatational wave speed field to auxiliary fields.
    void addVp(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AbsorbingDampersAuxiliaryFactory(const AbsorbingDampersAuxiliaryFactory &); ///< Not implemented.
    const AbsorbingDampersAuxiliaryFactory& operator=(const AbsorbingDampersAuxiliaryFactory&); ///< Not implemented

}; // class AbsorbingDampersAuxiliaryFactory

#endif // pylith_bc_absorbingdampersauxiliaryfactory_hh

// End of file
