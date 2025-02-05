// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "GravityRefState2D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

// ------------------------------------------------------------------------------------------------
namespace pylith {
    class _GravityRefState2D;
}
class pylith::_GravityRefState2D {
private:

    static const double GACC;
    static const double YMAX;

    // Density
    static double density(const double x,
                          const double y) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y) {
        return 3000.0;
    } // vs

    static const char* vs_units(void) {
        return "m/s";
    } // vs_units

    // Vp
    static double vp(const double x,
                     const double y) {
        return sqrt(3.0)*vs(x,y);
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    static double referenceMeanStress(const double x,
                                      const double y) {
        return density(x,y) * GACC * (y-YMAX);
    } // referenceMeanStress

    static double referenceShearStress(const double x,
                                       const double y) {
        return 0.0;
    } // referenceShearStress

    static const char* stress_units(void) {
        return "Pa";
    } // stress_units

    static double referenceStrain(const double x,
                                  const double y) {
        return 0.0;
    } // referencStrain

    static const char* strain_units(void) {
        return "none";
    } // strain_units

    static double setGravityAcc_x(const double x,
                                  const double y) {
        return 0.0;
    } // setGravityAcc_x

    static double setGravityAcc_y(const double x,
                                  const double y) {
        return -GACC;
    } // setGravityAcc_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units

    // Solution subfields (nondimensional).

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0;
    } // disp_y

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(s);
        assert(x);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

public:

    static
    TestLinearElasticity_Data* createData(void) {
        TestLinearElasticity_Data* data = new TestLinearElasticity_Data();assert(data);

        data->journalName = "GravityRefState2D";
        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(1.0e+03);
        data->normalizer.setTimeScale(2.0);
        data->normalizer.setPressureScale(2.25e+10);
        data->normalizer.computeDensityScale();

        delete data->gravityField;data->gravityField = new spatialdata::spatialdb::GravityField();
        data->gravityField->setGravityDir(0.0, -1.0, 0.0);
        data->gravityField->setGravityAcc(GACC);

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 6;
        static const char* _auxSubfields[6] = { // order must match order of subfields in auxiliary field
            "density",
            "gravitational_acceleration",
            "reference_stress",
            "reference_strain",
            "shear_modulus",
            "bulk_modulus",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(1, 1), // reference_stress
            pylith::topology::Field::Discretization(0, 1), // reference_strain
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.addValue("reference_stress_xx", referenceMeanStress, stress_units());
        data->auxDB.addValue("reference_stress_yy", referenceMeanStress, stress_units());
        data->auxDB.addValue("reference_stress_zz", referenceMeanStress, stress_units());
        data->auxDB.addValue("reference_stress_xy", referenceShearStress, stress_units());
        data->auxDB.addValue("reference_strain_xx", referenceStrain, strain_units());
        data->auxDB.addValue("reference_strain_yy", referenceStrain, strain_units());
        data->auxDB.addValue("reference_strain_zz", referenceStrain, strain_units());
        data->auxDB.addValue("reference_strain_xy", referenceStrain, strain_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->material.useBodyForce(false);
        data->rheology.useReferenceState(true);

        data->material.setIdentifier("elasticity");
        data->material.setName("material-id=24");
        data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        data->bcs.resize(1);
        pylith::bc::DirichletUserFn*bc = new pylith::bc::DirichletUserFn();assert(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        data->bcs[0] = bc;

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[1] = {
            solnkernel_disp,
        };
        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // _GravityRefState2D
const double pylith::_GravityRefState2D::GACC = 9.80665;
const double pylith::_GravityRefState2D::YMAX = +4.0e+3;

// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::GravityRefState2D::TriP1(void) {
    TestLinearElasticity_Data* data = pylith::_GravityRefState2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(1, 1), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
        pylith::topology::Field::Discretization(0, 1), // density
        pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
        pylith::topology::Field::Discretization(1, 1), // reference_stress
        pylith::topology::Field::Discretization(0, 1), // reference_strain
        pylith::topology::Field::Discretization(0, 1), // shear_modulus
        pylith::topology::Field::Discretization(0, 1), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::GravityRefState2D::TriP2(void) {
    TestLinearElasticity_Data* data = pylith::_GravityRefState2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
        pylith::topology::Field::Discretization(1, 2), // reference_stress
        pylith::topology::Field::Discretization(0, 2), // reference_strain
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::GravityRefState2D::TriP3(void) {
    TestLinearElasticity_Data* data = pylith::_GravityRefState2D::createData();assert(data);

    data->meshFilename = "data/tri.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
        pylith::topology::Field::Discretization(1, 3), // reference_stress
        pylith::topology::Field::Discretization(0, 3), // reference_strain
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TriP3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::GravityRefState2D::QuadQ1(void) {
    TestLinearElasticity_Data* data = pylith::_GravityRefState2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(1, 1), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
        pylith::topology::Field::Discretization(0, 1), // density
        pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
        pylith::topology::Field::Discretization(1, 1), // reference_stress
        pylith::topology::Field::Discretization(0, 1), // reference_strain
        pylith::topology::Field::Discretization(0, 1), // shear_modulus
        pylith::topology::Field::Discretization(0, 1), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::GravityRefState2D::QuadQ2(void) {
    TestLinearElasticity_Data* data = pylith::_GravityRefState2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
        pylith::topology::Field::Discretization(1, 2), // reference_stress
        pylith::topology::Field::Discretization(0, 2), // reference_strain
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::GravityRefState2D::QuadQ3(void) {
    TestLinearElasticity_Data* data = pylith::_GravityRefState2D::createData();assert(data);

    data->meshFilename = "data/quad.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[6] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
        pylith::topology::Field::Discretization(1, 3), // reference_stress
        pylith::topology::Field::Discretization(0, 3), // reference_strain
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // QuadQ3


// End of file
