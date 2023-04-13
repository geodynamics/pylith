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

#include <portinfo>

#include "AuxiliaryFactoryElastoplastic.hh" // implementation of object methods

#include "Material.hh"                            // USES Material
#include "Query.hh"                               // USES Query
#include "IsotropicDruckerPragerElastoplastic.hh" // USES IsotropicDruckerPragerElastoplastic

#include "pylith/topology/Field.hh"      // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh"   // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh"    // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryElastoplastic::AuxiliaryFactoryElastoplastic(void)
{
    GenericComponent::setName("auxiliaryfactoryelastoplastic");
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryElastoplastic::~AuxiliaryFactoryElastoplastic(void) {}

// ---------------------------------------------------------------------------------------------------------------------
// Add alpha yield subfield for Drucker-Prager elastoplastic model to auxiliary fields.
void pylith::materials::AuxiliaryFactoryElastoplastic::addAlphaYieldDruckerPrager(const PylithInt fitMohrCoulomb)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addAlphaYieldDruckerPrager(void)");

    const char *subfieldName = "alpha_yield";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    switch (fitMohrCoulomb)
    {
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_INSCRIBED:
    {
        pylith::materials::Query::alphaYieldDruckerPragerInscribedFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_INSCRIBED
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_MIDDLE:
    {
        pylith::materials::Query::alphaYieldDruckerPragerMiddleFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_MIDDLE
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_CIRCUMSCRIBED:
    {
        pylith::materials::Query::alphaYieldDruckerPragerCircumscribedFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_CIRCUMSCRIBED
    default:
        assert(0);
        throw std::logic_error("Unknown Mohr-Coulomb fit.");
        break;
    }

    PYLITH_METHOD_END;
} // addAlphaYieldDruckerPrager

// ---------------------------------------------------------------------------------------------------------------------
// Add alpha flow subfield for Drucker-Prager elastoplastic model to auxiliary fields.
void pylith::materials::AuxiliaryFactoryElastoplastic::addAlphaFlowDruckerPrager(const PylithInt fitMohrCoulomb)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addAlphaFlowDruckerPrager(void)");

    const char *subfieldName = "alpha_flow";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    switch (fitMohrCoulomb)
    {
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_INSCRIBED:
    {
        pylith::materials::Query::alphaFlowDruckerPragerInscribedFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_INSCRIBED
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_MIDDLE:
    {
        pylith::materials::Query::alphaFlowDruckerPragerMiddleFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_MIDDLE
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_CIRCUMSCRIBED:
    {
        pylith::materials::Query::alphaFlowDruckerPragerCircumscribedFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_CIRCUMSCRIBED
    default:
        assert(0);
        throw std::logic_error("Unknown Mohr-Coulomb fit.");
        break;
    }

    PYLITH_METHOD_END;
} // addAlphaFlowDruckerPrager

// ---------------------------------------------------------------------------------------------------------------------
// Add beta subfield for Drucker-Prager elastoplastic model to auxiliary fields.
void pylith::materials::AuxiliaryFactoryElastoplastic::addBetaDruckerPrager(const PylithInt fitMohrCoulomb)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBetaDruckerPrager(void)");

    const char *subfieldName = "beta";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    switch (fitMohrCoulomb)
    {
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_INSCRIBED:
    {
        pylith::materials::Query::betaDruckerPragerInscribedFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_INSCRIBED
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_MIDDLE:
    {
        pylith::materials::Query::betaDruckerPragerMiddleFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_MIDDLE
    case pylith::materials::IsotropicDruckerPragerElastoplastic::MOHR_COULOMB_CIRCUMSCRIBED:
    {
        pylith::materials::Query::betaDruckerPragerCircumscribedFromVM(subfieldName, this);
        break;
    } // MOHR_COULOMB_CIRCUMSCRIBED
    default:
        assert(0);
        throw std::logic_error("Unknown Mohr-Coulomb fit.");
        break;
    }

    PYLITH_METHOD_END;
} // addBetaDruckerPrager

// ---------------------------------------------------------------------------------------------------------------------
// Add plastic strain subfield to auxiliary fields.
void pylith::materials::AuxiliaryFactoryElastoplastic::addPlasticStrain(void)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPlasticStrain(void)");

    const char *subfieldName = "plastic_strain";
    const char *componentNames[6] = {
        "plastic_strain_xx",
        "plastic_strain_yy",
        "plastic_strain_zz",
        "plastic_strain_xy",
        "plastic_strain_yz",
        "plastic_strain_xz"};
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4
                                                                   : 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    description.hasHistory = true;
    description.historySize = 1;
    for (int i = 0; i < strainSize; ++i)
    {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPlasticStrain

// End of file
