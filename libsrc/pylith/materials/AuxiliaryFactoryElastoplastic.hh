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

/** @file libsrc/materials/AuxiliaryFactoryElastoplastic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for elastoplastic materials.
 */

#if !defined(pylith_materials_auxiliaryfactoryelastoplastic_hh)
#define pylith_materials_auxiliaryfactoryelastoplastic_hh

#include "materialsfwd.hh"                             // forward declarations
#include "pylith/materials/AuxiliaryFactoryElastic.hh" // ISA AuxiliaryFactoryElastic

/// @brief C++ helper class for setting up auxiliary fields for materials.
class pylith::materials::AuxiliaryFactoryElastoplastic : public pylith::materials::AuxiliaryFactoryElastic
{
    friend class TestAuxiliaryFactoryElastoplastic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    AuxiliaryFactoryElastoplastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryElastoplastic(void);

    /** Add alpha yield subfield for Drucker-Prager elastoplastic model to auxiliary subfields.
     *
     * @param[in] fitMohrCoulomb Type of fit to Mohr-Coulomb yield surface.
     */
    void addAlphaYieldDruckerPrager(const PylithInt fitMohrCoulomb);

    /** Add alpha flow subfield for Drucker-Prager elastoplastic model to auxiliary subfields.
     *
     * @param[in] fitMohrCoulomb Type of fit to Mohr-Coulomb yield surface.
     */
    void addAlphaFlowDruckerPrager(const PylithInt fitMohrCoulomb);

    /** Add beta subfield for Drucker-Prager elastoplastic model to auxiliary subfields.
     *
     * @param[in] fitMohrCoulomb Type of fit to Mohr-Coulomb yield surface.
     */
    void addBetaDruckerPrager(const PylithInt fitMohrCoulomb);

    /// Add plastic strain subfield to auxiliary subfields.
    void addPlasticStrain(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:
    AuxiliaryFactoryElastoplastic(const AuxiliaryFactoryElastoplastic &);                  ///< Not implemented.
    const AuxiliaryFactoryElastoplastic &operator=(const AuxiliaryFactoryElastoplastic &); ///< Not implemented

}; // class AuxiliaryFactoryElastoplastic

#endif // pylith_materials_auxiliaryfactoryelastoplastic_hh

// End of file
