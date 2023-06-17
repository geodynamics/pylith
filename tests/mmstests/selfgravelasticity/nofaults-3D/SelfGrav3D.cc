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

#include "SelfGrav3D.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh"         // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh"         // USES pythia::journal::debug_t

#include <cmath>

namespace pylith {
    class _SelfGrav3D;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::_SelfGrav3D {
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double BODYFORCE;
    static const double XMAX;

    static double sgn(const double x) {
        return signbit(x) < 0.0 ? -1.0 : 1.0;
    }

    static void xyzToRTP(double* r,
                         double* theta,
                         double* phi,
                         const double x,
                         const double y,
                         const double z) {
        *r = pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 0.5);
        *theta = *r > 0 ? acos(z / *r) : 0.0;
        *phi = x*x + y*y > 0.0 ? acos(x / sqrt(x * x + y * y)) * sgn(y) : 0.0;
    }

    static double gravConstN(void) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / velocityScale;
        const double gravScale = 1.0 / (densityScale * pow(TIMESCALE, 2));
        return 6.67e-11 / gravScale;
    } // Universal Gravitational Constant

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).

    // Density
    static double density(const double x,
                          const double y,
                          const double z) {
        return 2500.0;
    } // density

    static const char *density_units(void) {
        return "kg/m**3";
    } // density_units

    // Vs
    static double vs(const double x,
                     const double y,
                     const double z) {
        return 3000.0;
    } // vs

    static const char *vs_units(void) {
        return "m/s";
    } // vs_units

    // Vp
    static double vp(const double x,
                     const double y,
                     const double z) {
        return sqrt(3.0) * vs(x, y, z);
    } // vp

    static const char *vp_units(void) {
        return "m/s";
    } // vp_units

    static double bodyforce_x(const double x,
                              const double y,
                              const double z) {
        return BODYFORCE;
    } // bodyforce_x

    static double bodyforce_y(const double x,
                              const double y,
                              const double z) {
        return 0.0;
    } // bodyforce_y

    static double bodyforce_z(const double x,
                              const double y,
                              const double z) {
        return 0.0;
    } // bodyforce_z

    static double OuterRad(const double x,
                           const double y,
                           const double z) {
        return 1 * pow(10, 6);
    } // Radius in m (1000 km)

    static const char *bodyforce_units(void) {
        return "kg/(m**2*s**2)";
    } // bodyforce_units

    // Solution subfields (nondimensional)

    // Displacement

    static double disp_x(const double x,
                         const double y,
                         const double z) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / velocityScale;
        const double densityN = vs(x, y, z) / densityScale;
        const double vsN = vs(x, y, z) / velocityScale;
        const double vpN = vp(x, y, z) / velocityScale;
        const double outerRN = OuterRad(x, y, z) / LENGTHSCALE;

        // Put in spherical
        double r = 0.0;
        double theta = 0.0;
        double phi = 0.0;
        xyzToRTP(&r, &theta, &phi, x, y, z);

        const double fac1 = vsN * vsN / ((vpN + vsN) * (vpN - vsN));
        const double fac1_2 = vsN * vsN * densityN * (3 * vpN * vpN - 4 * vsN * vsN) / (vpN * vpN - vsN * vsN);
        const double fac2 = (5 * vpN * vpN - 4 * vsN * vsN) / (2 * (vpN * vpN - vsN * vsN));
        const double fac3 = (3 * vpN * vpN - 2 * vsN * vsN - 2 * vsN) / (2 * (vpN * vpN - vsN * vsN));

        const double ur = (2.0 / 15.0) * M_PI * gravConstN() * densityN * (fac1 / fac1_2) * (fac2 * outerRN * outerRN * r - fac3 * r * r * r);

        // put back in cartesian..

        const double x_disp = sin(theta) * cos(phi) * ur;

        return x_disp;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / velocityScale;
        const double densityN = vs(x, y, z) / densityScale;
        const double vsN = vs(x, y, z) / velocityScale;
        const double vpN = vp(x, y, z) / velocityScale;
        const double outerRN = OuterRad(x, y, z) / LENGTHSCALE;

        // Put in spherical
        double r = 0.0;
        double theta = 0.0;
        double phi = 0.0;
        xyzToRTP(&r, &theta, &phi, x, y, z);

        const double fac1 = vsN * vsN / ((vpN + vsN) * (vpN - vsN));
        const double fac1_2 = vsN * vsN * densityN * (3 * vpN * vpN - 4 * vsN * vsN) / (vpN * vpN - vsN * vsN);
        const double fac2 = (5 * vpN * vpN - 4 * vsN * vsN) / (2 * (vpN * vpN - vsN * vsN));
        const double fac3 = (3 * vpN * vpN - 2 * vsN * vsN - 2 * vsN) / (2 * (vpN * vpN - vsN * vsN));

        const double ur = (2.0 / 15.0) * M_PI * gravConstN() * densityN * (fac1 / fac1_2) * (fac2 * outerRN * outerRN * r - fac3 * r * r * r);
        // put back in cartesian..

        const double y_disp = sin(theta) * sin(phi) * ur;
        return y_disp;
    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / velocityScale;
        const double densityN = vs(x, y, z) / densityScale;
        const double vsN = vs(x, y, z) / velocityScale;
        const double vpN = vp(x, y, z) / velocityScale;
        const double outerRN = OuterRad(x, y, z) / LENGTHSCALE;

        // Put in spherical
        double r = 0.0;
        double theta = 0.0;
        double phi = 0.0;
        xyzToRTP(&r, &theta, &phi, x, y, z);

        const double fac1 = vsN * vsN / ((vpN + vsN) * (vpN - vsN));
        const double fac1_2 = vsN * vsN * densityN * (3 * vpN * vpN - 4 * vsN * vsN) / (vpN * vpN - vsN * vsN);
        const double fac2 = (5 * vpN * vpN - 4 * vsN * vsN) / (2 * (vpN * vpN - vsN * vsN));
        const double fac3 = (3 * vpN * vpN - 2 * vsN * vsN - 2 * vsN) / (2 * (vpN * vpN - vsN * vsN));

        const double ur = (2.0 / 15.0) * M_PI * gravConstN() * densityN * (fac1 / fac1_2) * (fac2 * outerRN * outerRN * r - fac3 * r * r * r);

        // put back in cartesian..
        const double z_disp = cos(theta) * ur;

        return z_disp;
    } // disp_z

    // potential
    static double potential(const double x,
                            const double y,
                            const double z) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double densityScale = PRESSURESCALE / velocityScale;
        const double densityN = vs(x, y, z) / densityScale;

        const double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

        return (4 / 3) * M_PI * gravConstN() * densityN * r;
    } // potential

    static PetscErrorCode solnkernel_potential(PetscInt spaceDim,
                                               PetscReal t,
                                               const PetscReal x[],
                                               PetscInt numComponents,
                                               PetscScalar *s,
                                               void *context) {
        assert(3 == spaceDim);
        assert(x);
        assert(1 == numComponents);
        assert(s);

        s[0] = potential(x[0], x[1], x[2]);

        return PETSC_SUCCESS;
    } // solnkernel_potential

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar *s,
                                          void *context) {
        assert(3 == spaceDim);
        assert(3 == numComponents);
        assert(s);

        s[0] = disp_x(x[0], x[1], x[2]);
        s[1] = disp_y(x[0], x[1], x[2]);
        s[2] = disp_z(x[0], x[1], x[2]);

        return 0;
    } // solnkernel_disp

public:

    static TestSelfGrav_Data *createData(void) {
        TestSelfGrav_Data *data = new TestSelfGrav_Data();
        assert(data);

        data->journalName = "SelfGrav3D";

        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(LENGTHSCALE);
        data->normalizer.setTimeScale(TIMESCALE);
        data->normalizer.setPressureScale(PRESSURESCALE);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 4;
        static const char *_auxSubfields[4] = {
            // order must match order of subfields in auxiliary field
            "density",
            "body_force",
            "shear_modulus",
            "bulk_modulus",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // body_force
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

        data->auxDB.addValue("density", density, density_units());
        data->auxDB.addValue("vp", vp, vp_units());
        data->auxDB.addValue("vs", vs, vs_units());
        data->auxDB.addValue("body_force_x", bodyforce_x, bodyforce_units());
        data->auxDB.addValue("body_force_y", bodyforce_y, bodyforce_units());
        data->auxDB.addValue("body_force_z", bodyforce_z, bodyforce_units());
        data->auxDB.setCoordSys(data->cs);

        data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
        data->material.useBodyForce(true);
        data->rheology.useReferenceState(false);

        data->material.setDescription("Isotropic Linear Elascitity Plane Strain");
        data->material.setLabelValue(24);

        static const PylithInt constrainedDOF[3] = {0, 1, 2};
        static const PylithInt numConstrained = 3;
        data->bcs.resize(1);
        pylith::bc::DirichletUserFn *bc = new pylith::bc::DirichletUserFn();
        assert(bc);
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(constrainedDOF, numConstrained);
        bc->setUserFn(solnkernel_disp);
        data->bcs[0] = bc;

        { // potential
            static const PylithInt constrainedDOF[1] = {0};
            static const PylithInt numConstrainedDOF = 1;
            pylith::bc::DirichletUserFn *bc = new pylith::bc::DirichletUserFn();
            assert(bc);
            bc->setSubfieldName("potential");
            bc->setLabelName("boundary");
            bc->setLabelValue(1);
            bc->setConstrainedDOF(constrainedDOF, numConstrainedDOF);
            bc->setUserFn(solnkernel_potential);
            data->bcs[1] = bc;
        } // potential

        static const pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
            solnkernel_disp,
            solnkernel_potential,
        };

        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn *>(_exactSolnFns);
        data->exactSolnDotFns = nullptr;

        return data;
    } // createData

}; // SelfGrav3D
const double pylith::_SelfGrav3D::LENGTHSCALE = 1.0e+3;
const double pylith::_SelfGrav3D::TIMESCALE = 2.0;
const double pylith::_SelfGrav3D::PRESSURESCALE = 2.25e+10;
const double pylith::_SelfGrav3D::BODYFORCE = 5.0e+3;
const double pylith::_SelfGrav3D::XMAX = 4.0e+3;

// ------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------
pylith::TestSelfGrav_Data *
pylith::SelfGrav3D::TetP2(void) {
    TestSelfGrav_Data *data = pylith::_SelfGrav3D::createData();
    assert(data);

    data->meshFilename = "data/tet.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(2, 2), // potential
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // body_force
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization *>(_auxDiscretizations);

    return data;
} // TetP3


// End of file
