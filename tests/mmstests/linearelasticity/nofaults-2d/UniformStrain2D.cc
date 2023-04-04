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

#include "TestLinearElasticity.hh" // ISA TestIsotropicLinearElasticity2D

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

namespace pylith {
    namespace mmstests {
        class UniformStrain2D;

        class UniformStrain2D_TriP1;
        class UniformStrain2D_TriP2;
        class UniformStrain2D_TriP3;
        class UniformStrain2D_TriP4;

        class UniformStrain2D_QuadQ1;
        class UniformStrain2D_QuadQ1Distorted;
        class UniformStrain2D_QuadQ2;
        class UniformStrain2D_QuadQ3;
        class UniformStrain2D_QuadQ4;
    } // tests/mmstests
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D :
    public pylith::mmstests::TestLinearElasticity {
    // Spatial database user functions for auxiiliary subfields (includes derived fields).

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

    // Solution subfields.

    static double strain_xx(void) {
        return 0.1;
    } // strain_xx

    static double strain_yy(void) {
        return 0.25;
    } // strain_yy

    static double strain_xy(void) {
        return 0.3;
    } // strain_xy

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return strain_xx()*x + strain_xy()*y;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return strain_xy()*x + strain_yy()*y;
    } // disp_y

    static const char* disp_units(void) {
        return "m";
    } // disp_units

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

protected:

    void setUp(void) {
        TestLinearElasticity::setUp();

        // Overwrite component names for control of journals at test level.
        GenericComponent::setName("UniformStrain2D");

        CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        _data->normalizer.setLengthScale(1.0e+03);
        _data->normalizer.setTimeScale(2.0);
        _data->normalizer.setPressureScale(2.25e+10);
        _data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        _data->numAuxSubfields = 3;
        static const char* _auxSubfields[3] = {"density", "shear_modulus", "bulk_modulus"};
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization const*>(_auxDiscretizations);

        _data->auxDB.addValue("density", density, density_units());
        _data->auxDB.addValue("vp", vp, vp_units());
        _data->auxDB.addValue("vs", vs, vs_units());
        _data->auxDB.setCoordSys(_data->cs);

        _data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        _data->material.useBodyForce(false);
        _data->rheology.useReferenceState(false);

        _data->material.setDescription("Isotropic Linear Elasticity Plane Strain");
        _data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[2] = {0, 1};
        static const PylithInt numConstrained = 2;
        _data->bc.setSubfieldName("displacement");
        _data->bc.setLabelName("boundary");
        _data->bc.setLabelValue(1);
        _data->bc.setConstrainedDOF(constrainedDOF, numConstrained);
        _data->bc.setUserFn(solnkernel_disp);

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        const pylith::topology::Field* solution = _problem->getSolution();
        CPPUNIT_ASSERT(solution);

        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(solution->getDM(), &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // UniformStrain2D

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_TriP1 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_TriP1,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.msh";
        _data->useAsciiMesh = false;

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_TriP2 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_TriP2,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_TriP2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_TriP3 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_TriP3,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.msh";
        _data->useAsciiMesh = false;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_TriP3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_TriP4 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_TriP4,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_TriP4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_TriP4);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_QuadQ1 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_QuadQ1,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_QuadQ1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_QuadQ1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_QuadQ1Distorted :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_QuadQ1Distorted,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_distorted.mesh";

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(1, 1), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_QuadQ1Distorted
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_QuadQ1Distorted);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_QuadQ2 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_QuadQ2,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.msh";
        _data->useAsciiMesh = false;

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(2, 2), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_QuadQ2
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_QuadQ3 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_QuadQ3,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(3, 3), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_QuadQ3
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_QuadQ3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::UniformStrain2D_QuadQ4 :
    public pylith::mmstests::UniformStrain2D {
    CPPUNIT_TEST_SUB_SUITE(UniformStrain2D_QuadQ4,
                           TestLinearElasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        UniformStrain2D::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad.mesh";

        static const pylith::topology::Field::Discretization _auxDiscretizations[3] = {
            pylith::topology::Field::Discretization(0, 4), // density
            pylith::topology::Field::Discretization(0, 4), // shear_modulus
            pylith::topology::Field::Discretization(0, 4), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        _data->numSolnSubfields = 1;
        static const pylith::topology::Field::Discretization _solnDiscretizations[1] = {
            pylith::topology::Field::Discretization(4, 4), // disp
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    } // setUp

}; // UniformStrain2D_QuadQ4
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::UniformStrain2D_QuadQ4);

// End of file
