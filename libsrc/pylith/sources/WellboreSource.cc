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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/sources/WellboreSource.hh" // implementation of object methods

#include "pylith/sources/AuxiliaryFactoryWellboreSource.hh" // USES AuxiliaryFactoryWellboreSource
#include "pylith/feassemble/IntegratorDomain.hh"            // USES IntegratorDomain
#include "pylith/topology/Mesh.hh"                          // USES Mesh
#include "pylith/topology/Field.hh"                         // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh"                      // USES FieldOps

#include "pylith/fekernels/WellboreSource.hh" // USES WellboreSource kernels

#include "pylith/utils/error.hh"    // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/geocoords/CoordSys.hh"   // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::WellboreSource::WellboreSource(void) : _auxiliaryFactory(new pylith::sources::AuxiliaryFactoryWellboreSource)
{
    pylith::utils::PyreComponent::setName("wellboresource");
} // constructor

// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::WellboreSource::~WellboreSource(void)
{
    deallocate();
} // destructor

// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::sources::WellboreSource::deallocate(void)
{
    Source::deallocate();

    delete _auxiliaryFactory;
    _auxiliaryFactory = NULL;
} // deallocate

// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void pylith::sources::WellboreSource::verifyConfiguration(const pylith::topology::Field &solution) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution=" << solution.getLabel() << ")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("pressure"))
    {
        throw std::runtime_error("Cannot find 'pressure' field in solution; required for 'WellboreSource'.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator *
pylith::sources::WellboreSource::createIntegrator(const pylith::topology::Field &solution)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution=" << solution.getLabel() << ")");

    printf("In WellboreSource begin\n");
    PetscErrorCode err;
    PetscDM dmSoln = solution.getDM();
    assert(dmSoln);
    // transform points of source to mesh coordinates in python
    // DM from solution
    Vec vecPoints;
    DMLabel label;
    PetscSF sfPoints = NULL;
    const PetscInt *localPoints;
    const PetscSFNode *remotePoints;
    PetscInt numRoots = -1, numLeaves, dim, d;
    PetscMPIInt rank;
    PetscScalar *a;

    err = DMGetCoordinateDim(dmSoln, &dim);
    PYLITH_CHECK_ERROR(err);
    err = VecCreateMPIWithArray(PetscObjectComm((PetscObject)dmSoln), dim, _pointCoords.size(), PETSC_DECIDE,
                                &_pointCoords[0], &vecPoints);
    PYLITH_CHECK_ERROR(err);

    // Debug
    PetscPrintf(PetscObjectComm((PetscObject)dmSoln), "_pointCoords\n");
    PetscPrintf(PetscObjectComm((PetscObject)dmSoln), " x = %g\n", _pointCoords[0]);
    PetscPrintf(PetscObjectComm((PetscObject)dmSoln), " y = %g\n", _pointCoords[1]);
    // Erzatz from ex17
    // err = VecCreateSeq(PETSC_COMM_SELF, dim, &vecPoints);PYLITH_CHECK_ERROR(err);
    // err = VecSetBlockSize(vecPoints, _pointCoords.size());PYLITH_CHECK_ERROR(err);
    // err = VecGetArray(vecPoints, &a);PYLITH_CHECK_ERROR(err);
    // for (d = 0; d < _pointCoords.size(); ++d) {
    //     a[d] = _pointCoords[d];
    // }
    // err = VecRestoreArray(vecPoints, &a);PYLITH_CHECK_ERROR(err);

    err = DMLocatePoints(dmSoln, vecPoints, DM_POINTLOCATION_NONE, &sfPoints);
    PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&vecPoints);
    PYLITH_CHECK_ERROR(err);
    err = DMCreateLabel(dmSoln, PyreComponent::getIdentifier());
    PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmSoln, PyreComponent::getIdentifier(), &label);
    PYLITH_CHECK_ERROR(err);
    err = PetscSFGetGraph(sfPoints, &numRoots, &numLeaves, &localPoints, &remotePoints);
    PYLITH_CHECK_ERROR(err);
    err = MPI_Comm_rank(PetscObjectComm((PetscObject)dmSoln), &rank);
    PYLITH_CHECK_ERROR(err);
    // Debug
    PetscPrintf(PetscObjectComm((PetscObject)dmSoln), "localPoints: %D\n", numLeaves);
    for (PetscInt p = 0; p < numLeaves; ++p)
    {
        if (remotePoints[p].rank == rank)
        {
            err = DMLabelSetValue(label, remotePoints[p].index, 2);
            PYLITH_CHECK_ERROR(err);
        }
    } // for
    err = PetscSFDestroy(&sfPoints);
    PYLITH_CHECK_ERROR(err);

    pylith::feassemble::IntegratorDomain *integrator = new pylith::feassemble::IntegratorDomain(this);
    assert(integrator);
    integrator->setLabelName(PyreComponent::getIdentifier());
    integrator->setLabelValue(_sourceId);
    printf("In WellboreSource end\n");
    DMView(dmSoln, NULL);

    _setKernelsResidual(integrator, solution);
    _setKernelsJacobian(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator

// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field *
pylith::sources::WellboreSource::createAuxiliaryField(const pylith::topology::Field &solution,
                                                      const pylith::topology::Mesh &domainMesh)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution=" << solution.getLabel() << ", domainMesh=" << typeid(domainMesh).name() << ")");

    pylith::topology::Field *auxiliaryField = new pylith::topology::Field(domainMesh);
    assert(auxiliaryField);
    auxiliaryField->setLabel("WellboreSource auxiliary field");

    assert(_normalizer);
    _auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.getDimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    // add in aux specific to peaceman
    _auxiliaryFactory->addFluidDensity();          // 0
    _auxiliaryFactory->addFluidViscosity();        // 1
    _auxiliaryFactory->addIsotropicPermeability(); // 2
    _auxiliaryFactory->addWellboreRadius();        // 3
    _auxiliaryFactory->addWellboreLength();        // 4
    _auxiliaryFactory->addWellborePressure();      // 5
    _auxiliaryFactory->addWellboreCharacter();     // 6
    _auxiliaryFactory->addElementDimensions();     // 7

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setValuesFromDB();

    // Debug option
    auxiliaryField->view("Wellbore auxiliary field.");

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField

// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field *
pylith::sources::WellboreSource::createDerivedField(const pylith::topology::Field &solution,
                                                    const pylith::topology::Mesh &domainMesh)
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution=" << solution.getLabel() << ", domainMesh=)" << typeid(domainMesh).name() << ") empty method");

    PYLITH_METHOD_RETURN(NULL);
} // createDerivedField

// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory *
pylith::sources::WellboreSource::_getAuxiliaryFactory(void)
{
    return _auxiliaryFactory;
} // _getAuxiliaryFactory

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for residual.
void pylith::sources::WellboreSource::_setKernelsResidual(pylith::feassemble::IntegratorDomain *integrator,
                                                          const topology::Field &solution) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator=" << integrator << ", solution=" << solution.getLabel() << ")");

    const spatialdata::geocoords::CoordSys *coordsys = solution.getMesh().getCoordSys();

    std::vector<ResidualKernels> kernels;

    switch (_formulation)
    {
    case QUASISTATIC:
    {
        // Pressure
        const PetscPointFunc f0p = pylith::fekernels::WellboreSource::f0p;
        const PetscPointFunc f1p = NULL;

        kernels.resize(1);
        kernels[0] = ResidualKernels("pressure", pylith::feassemble::Integrator::RESIDUAL_LHS, f0p, f1p);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX:
    {
        break;
    } // DYNAMIC
    case DYNAMIC:
    {
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void pylith::sources::WellboreSource::_setKernelsJacobian(pylith::feassemble::IntegratorDomain *integrator,
                                                          const topology::Field &solution) const
{
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator=" << integrator << ", solution=" << solution.getLabel() << ")");

    const spatialdata::geocoords::CoordSys *coordsys = solution.getMesh().getCoordSys();

    std::vector<JacobianKernels> kernels;

    switch (_formulation)
    {
    case QUASISTATIC:
    {
        const PetscPointJac Jf0pp = pylith::fekernels::WellboreSource::Jf0pp;
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = NULL;

        kernels.resize(1);
        kernels[0] = JacobianKernels("pressure", "pressure", pylith::feassemble::Integrator::JACOBIAN_LHS, Jf0pp, Jf1pp, Jf2pp, Jf3pp);
        break;
    } // QUASISTATIC
    case DYNAMIC:
    case DYNAMIC_IMEX:
    {
        break;
    } // DYNAMIC_IMEX
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    assert(integrator);
    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian

// End of file
