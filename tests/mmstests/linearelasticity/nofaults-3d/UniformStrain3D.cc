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

#include "UniformStrain3D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

namespace pylith {
    class _UniformStrain3D;
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::_UniformStrain3D {
private:

    // Density
    static double density(const double x,
                          const double y,
                          const double z) {
        return 2500.0;
    } // density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y,
                     const double z) {
        return 3000.0;
    } // vs

    static const char* vs_units(void) {
        return "m/s";
    } // vs_units

    // Vp
    static double vp(const double x,
                     const double y,
                     const double z) {
        return sqrt(3.0)*vs(x,y,z);
    } // vp

    static const char* vp_units(void) {
        return "m/s";
    } // vp_units

    // Solution subfields.

    static double strain_xx(void) {
        return 0.1;
    } // strain_xx

    static double strain_yy(void) {
        return 0.25;
    } // strain_yy

    static double strain_zz(void) {
        return -0.05;
    } // strain_zz

    static double strain_xy(void) {
        return 0.3;
    } // strain_xy

    static double strain_yz(void) {
        return 0;
    } // strain_yz

    static double strain_xz(void) {
        return -0.08;
    } // strain_xz

    // Displacement
    static double disp_x(const double x,
                         const double y,
                         const double z) {
        return strain_xx()*x + strain_xy()*y + strain_xz()*z;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {
        return strain_xy()*x + strain_yy()*y * strain_yz()*z;
    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {
        return strain_xz()*x + strain_yz()*y * strain_zz()*z;
    } // disp_z

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        assert(3 == spaceDim);
        assert(x);
        assert(3 == numComponents);
        assert(s);

        s[0] = disp_x(x[0], x[1], x[2]);
        s[1] = disp_y(x[0], x[1], x[2]);
        s[2] = disp_z(x[0], x[1], x[2]);

        return 0;
    } // solnkernel_disp

public:

    static
    TestLinearElasticity_Data* createData(void) {
        TestLinearElasticity_Data* data = new TestLinearElasticity_Data();assert(data);

        data->journalName = "UniformStrain3D";

        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(1.0e+03);
        data->normalizer.setTimeScale(2.0);
        data->normalizer.setPressureScale(2.25e+10);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->material.useBodyForce(false);
        data->rheology.useReferenceState(false);

        data->material.setDescription("Isotropic Linear Elascitity");
        data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[3] = { 0, 1, 2 };
        static const PylithInt numConstrained = 3;
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

}; // UniformStrain3D

// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::TetP1(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/tet.mesh";

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(1, 1), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TetP1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::TetP2(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/tet.msh";
    data->useAsciiMesh = false;

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TetP2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::TetP3(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/tet.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TetP3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::TetP4(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/tet.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 4), // density
        pylith::topology::Field::Discretization(0, 4), // shear_modulus
        pylith::topology::Field::Discretization(0, 4), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(4, 4), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // TetP4


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::HexQ1(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/hex.msh";
    data->useAsciiMesh = false;

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(1, 1), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // HexQ1


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::HexQ2(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/hex.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(2, 2), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // HexQ2


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::HexQ3(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/hex.msh";
    data->useAsciiMesh = false;

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(3, 3), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // HexQ3


// ------------------------------------------------------------------------------------------------
pylith::TestLinearElasticity_Data*
pylith::UniformStrain3D::HexQ4(void) {
    TestLinearElasticity_Data* data = pylith::_UniformStrain3D::createData();assert(data);

    data->meshFilename = "data/hex.mesh";

    static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
        pylith::topology::Field::Discretization(0, 4), // density
        pylith::topology::Field::Discretization(0, 4), // shear_modulus
        pylith::topology::Field::Discretization(0, 4), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    data->numSolnSubfields = 1;
    static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
        pylith::topology::Field::Discretization(4, 4), // disp
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
} // HexQ4


// End of file
