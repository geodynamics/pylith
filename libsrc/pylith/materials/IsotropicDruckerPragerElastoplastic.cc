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

#include "pylith/materials/IsotropicDruckerPragerElastoplastic.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastoplastic.hh"       // USES AuxiliaryFactoryElastoplastic
#include "pylith/fekernels/IsotropicDruckerPragerElastoplastic.hh" // USES IsotropicDruckerPragerElastoplastic kernels
#include "pylith/feassemble/Integrator.hh"                         // USES Integrator
#include "pylith/utils/journals.hh"                                // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh"                                   // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicDruckerPragerElastoplastic::IsotropicDruckerPragerElastoplastic(void) : _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
                                                                                                    _useReferenceState(false),
                                                                                                    _fitMohrCoulomb(MOHR_COULOMB_INSCRIBED),
                                                                                                    _allowTensileYield(false)
{
    _lhsJacobianTriggers = pylith::feassemble::Integrator::NEW_JACOBIAN_ALWAYS;
    pylith::utils::PyreComponent::setName("isotropicdruckerpragerep");
} // constructor

// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicDruckerPragerElastoplastic::~IsotropicDruckerPragerElastoplastic(void)
{
    deallocate();
} // destructor

// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::materials::IsotropicDruckerPragerElastoplastic::deallocate(void)
{
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;
    _auxiliaryFactory = NULL;
} // deallocate

// ---------------------------------------------------------------------------------------------------------------------
// Set fit to Mohr-Coulomb surface.
void pylith::materials::IsotropicDruckerPragerElastoplastic::fitMohrCoulomb(const FitMohrCoulombEnum value)
{
    PYLITH_COMPONENT_DEBUG("fitMohrCoulomb=" << value << ")");

    _fitMohrCoulomb = value;
} // fitMohrCoulomb

// ---------------------------------------------------------------------------------------------------------------------
// Get fit to Mohr-Coulomb surface.
int pylith::materials::IsotropicDruckerPragerElastoplastic::fitMohrCoulomb(void) const
{
    return _fitMohrCoulomb;
} // fitMohrCoulomb

// ---------------------------------------------------------------------------------------------------------------------
// Set flag for whether to allow tensile yield.
void pylith::materials::IsotropicDruckerPragerElastoplastic::allowTensileYield(const bool value)
{
    PYLITH_COMPONENT_DEBUG("allowTensileYield=" << value << ")");

    _allowTensileYield = value;
} // allowTensileYield

// ---------------------------------------------------------------------------------------------------------------------
// Get flag for whether to allow tensile yield.
bool pylith::materials::IsotropicDruckerPragerElastoplastic::allowTensileYield(void) const
{
    return _allowTensileYield;
} // allowTensileYield

// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void pylith::materials::IsotropicDruckerPragerElastoplastic::useReferenceState(const bool value)
{
    PYLITH_COMPONENT_DEBUG("useReferenceState=" << value << ")");

    _useReferenceState = value;
} // useReferenceState

// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool pylith::materials::IsotropicDruckerPragerElastoplastic::useReferenceState(void) const
{
    return _useReferenceState;
} // useReferenceState

// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity *
pylith::materials::IsotropicDruckerPragerElastoplastic::getAuxiliaryFactory(void)
{
    return _auxiliaryFactory;
} // getAuxiliaryFactory

// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void pylith::materials::IsotropicDruckerPragerElastoplastic::addAuxiliarySubfields(void)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState)
    {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();
    _auxiliaryFactory->addAlphaYield(_fitMohrCoulomb);
    _auxiliaryFactory->addBeta(_fitMohrCoulomb);
    _auxiliaryFactory->addAlphaFlow(_fitMohrCoulomb);
    _auxiliaryFactory->addPlasticStrain();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicDruckerPragerElastoplastic::getKernelf1v(const spatialdata::geocoords::CoordSys *coordsys) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1v(coordsys=" << typeid(coordsys).name() << ")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::f1v_infinitesimalStrain_tensile : (!_useReferenceState && 2 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::f1v_infinitesimalStrain_tensile
                                                                                                                                                               : (_useReferenceState && 3 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::f1v_infinitesimalStrain_refState_tensile
                                                                                                                                                               : (_useReferenceState && 2 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::f1v_infinitesimalStrain_refState_tensile
                                                                                                                                                               : (!_useReferenceState && 3 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::f1v_infinitesimalStrain_noTensile
                                                                                                                                                               : (!_useReferenceState && 2 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::f1v_infinitesimalStrain_noTensile
                                                                                                                                                               : (_useReferenceState && 3 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::f1v_infinitesimalStrain_refState_noTensile
                                                                                                                                                               : (_useReferenceState && 2 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::f1v_infinitesimalStrain_refState_noTensile
                                                                                                                                                                                                                                : NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1v

// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicDruckerPragerElastoplastic::getKernelJf3vu(const spatialdata::geocoords::CoordSys *coordsys) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3vu(coordsys=" << typeid(coordsys).name() << ")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (!_useReferenceState && 3 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::Jf3vu_infinitesimalStrain_tensile : (!_useReferenceState && 2 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::Jf3vu_infinitesimalStrain_tensile
                                                                                                                                                                 : (_useReferenceState && 3 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::Jf3vu_infinitesimalStrain_refState_tensile
                                                                                                                                                                 : (_useReferenceState && 2 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::Jf3vu_infinitesimalStrain_refState_tensile
                                                                                                                                                                 : (!_useReferenceState && 3 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::Jf3vu_infinitesimalStrain_noTensile
                                                                                                                                                                 : (!_useReferenceState && 2 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::Jf3vu_infinitesimalStrain_noTensile
                                                                                                                                                                 : (_useReferenceState && 3 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::Jf3vu_infinitesimalStrain_refState_noTensile
                                                                                                                                                                 : (_useReferenceState && 2 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::Jf3vu_infinitesimalStrain_refState_noTensile
                                                                                                                                                                                                                                  : NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJf3vu

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicDruckerPragerElastoplastic::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys *coordsys) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys=" << typeid(coordsys).name() << ")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::cauchyStress_infinitesimalStrain_tensile_asVector : (!_useReferenceState && 2 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::cauchyStress_infinitesimalStrain_tensile_asVector
                                                                                                                                                                                 : (_useReferenceState && 3 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::cauchyStress_infinitesimalStrain_refState_tensile_asVector
                                                                                                                                                                                 : (_useReferenceState && 2 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::cauchyStress_infinitesimalStrain_refState_tensile_asVector
                                                                                                                                                                                 : (!_useReferenceState && 3 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::cauchyStress_infinitesimalStrain_noTensile_asVector
                                                                                                                                                                                 : (!_useReferenceState && 2 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::cauchyStress_infinitesimalStrain_noTensile_asVector
                                                                                                                                                                                 : (_useReferenceState && 3 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::cauchyStress_infinitesimalStrain_refState_noTensile_asVector
                                                                                                                                                                                 : (_useReferenceState && 2 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::cauchyStress_infinitesimalStrain_refState_noTensile_asVector
                                                                                                                                                                                                                                                  : NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector

// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void pylith::materials::IsotropicDruckerPragerElastoplastic::updateKernelConstants(pylith::real_array *kernelConstants,
                                                                                   const PylithReal dt) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants" << kernelConstants << ", dt=" << dt << ")");

    assert(kernelConstants);

    if (1 != kernelConstants->size())
    {
        kernelConstants->resize(1);
    }
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants

// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void pylith::materials::IsotropicDruckerPragerElastoplastic::addKernelsUpdateStateVars(std::vector<ProjectKernels> *kernels,
                                                                                       const spatialdata::geocoords::CoordSys *coordsys) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels=" << kernels << ", coordsys=" << coordsys << ")");

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc funcPlasticStrain =
        (!_useReferenceState && 3 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::plasticStrain_infinitesimalStrain_tensile_asVector : (!_useReferenceState && 2 == spaceDim && _allowTensileYield) ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::plasticStrain_infinitesimalStrain_tensile_asVector
                                                                                                                                                                                  : (_useReferenceState && 3 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::plasticStrain_infinitesimalStrain_refState_tensile_asVector
                                                                                                                                                                                  : (_useReferenceState && 2 == spaceDim && _allowTensileYield)    ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::plasticStrain_infinitesimalStrain_refState_tensile_asVector
                                                                                                                                                                                  : (!_useReferenceState && 3 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::plasticStrain_infinitesimalStrain_noTensile_asVector
                                                                                                                                                                                  : (!_useReferenceState && 2 == spaceDim && !_allowTensileYield)  ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::plasticStrain_infinitesimalStrain_noTensile_asVector
                                                                                                                                                                                  : (_useReferenceState && 3 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplastic3D::plasticStrain_infinitesimalStrain_refState_noTensile_asVector
                                                                                                                                                                                  : (_useReferenceState && 2 == spaceDim && !_allowTensileYield)   ? pylith::fekernels::IsotropicDruckerPragerElastoplasticPlaneStrain::plasticStrain_infinitesimalStrain_refState_noTensile_asVector
                                                                                                                                                                                                                                                   : NULL;
    assert(kernels);
    size_t prevNumKernels = kernels->size();
    kernels->resize(prevNumKernels + 1);
    (*kernels)[prevNumKernels + 0] = ProjectKernels("plastic_strain", funcPlasticStrain);

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars

// End of file
