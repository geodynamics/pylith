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
#include "math.h" // Implementation of math

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES pythia::journal::debug_t

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

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).

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

    static double GravConst(const double x,
                              const double y,
                              const double z) {
        return 6.67*pow(10,-11);
    } // Universal Gravitational Constant

    static double OuterRad(const double x,
                              const double y,
                              const double z) {
        return 1*pow(10,6);
    } // Radius in m (1000 km)



    static const char* bodyforce_units(void) {
        return "kg/(m**2*s**2)";
    } // bodyforce_units

    // Solution subfields (nondimensional)

    // Displacement

    // Define a function for signum 



    static double disp_x(const double x,
                         const double y,
                         const double z) {

        // Put in spherical

        const double r = pow(pow(x,2)+pow(x,2)+pow(x,2),0.5);
        const double fac1 = vs(x,y,z)*vs(x,y,z)/((vp(x,y,z)+vs(x,y,z))*(vp(x,y,z)-vs(x,y,z)))
        const double fac1_2 = vs(x,y,z)*vs(x,y,z)*density(x.y.z)*(3*vp(x,y,z)*vp(x,y,z)-4*vs(x,y,z)*vs(x,y,z))/(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z))
        const double fac2 = (5*vp(x,y,z)*vp(x,y,z)-4*vs(x,y,z)*vs(x,y,z))/(2*(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z)))
        const double fac3 = (3*vp(x,y,z)*vp(x,y,z)- 2*vs(x,y,z)*vs(x,y,z)-2*vs(x,y,z))/(2*(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z)))

        const double ur = (2/15)*3.14159*GravConst(x,y,z)*density(x,y,z)*(fac1/fac1_2)*(fac2*OuterRad(x,y,z)*OuterRad(x,y,z)*r-fac3*r*r*r)

        //put back in cartesian..


        int sgn(double d){ 
            return d<-eps?-1:d>eps;
        }
        
        
        theta = acos(z/(r))
        phi= acos(x/(x*x+y*y))*sgn(y)
        x_disp = sin(theta)*cos(phi)*ur
    


        return x_disp;
    } // disp_x

    static double disp_y(const double x,
                         const double y,
                         const double z) {

        // Put in spherical

        const double r = pow(pow(x,2)+pow(x,2)+pow(x,2),0.5);
        const double fac1 = vs(x,y,z)*vs(x,y,z)/((vp(x,y,z)+vs(x,y,z))*(vp(x,y,z)-vs(x,y,z)))
        const double fac1_2 = vs(x,y,z)*vs(x,y,z)*density(x.y.z)*(3*vp(x,y,z)*vp(x,y,z)-4*vs(x,y,z)*vs(x,y,z))/(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z))
        const double fac2 = (5*vp(x,y,z)*vp(x,y,z)-4*vs(x,y,z)*vs(x,y,z))/(2*(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z)))
        const double fac3 = (3*vp(x,y,z)*vp(x,y,z)- 2*vs(x,y,z)*vs(x,y,z)-2*vs(x,y,z))/(2*(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z)))

        const double ur = (2/15)*3.14159*GravConst(x,y,z)*density(x,y,z)*(fac1/fac1_2)*(fac2*OuterRad(x,y,z)*OuterRad(x,y,z)*r-fac3*r*r*r)

        //put back in cartesian..


        int sgn(double d){ 
            return d<-eps?-1:d>eps;
        }
        
        
        theta = acos(z/(r))
        phi= acos(x/(x*x+y*y))*sgn(y)
        y_disp = sin(theta)*sin(phi)*ur
        return y_disp;
    } // disp_y

    static double disp_z(const double x,
                         const double y,
                         const double z) {

        // Put in spherical


        const double r = pow(pow(x,2)+pow(x,2)+pow(x,2),0.5);
        const double fac1 = vs(x,y,z)*vs(x,y,z)/((vp(x,y,z)+vs(x,y,z))*(vp(x,y,z)-vs(x,y,z)))
        const double fac1_2 = vs(x,y,z)*vs(x,y,z)*density(x.y.z)*(3*vp(x,y,z)*vp(x,y,z)-4*vs(x,y,z)*vs(x,y,z))/(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z))
        const double fac2 = (5*vp(x,y,z)*vp(x,y,z)-4*vs(x,y,z)*vs(x,y,z))/(2*(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z)))
        const double fac3 = (3*vp(x,y,z)*vp(x,y,z)- 2*vs(x,y,z)*vs(x,y,z)-2*vs(x,y,z))/(2*(vp(x,y,z)*vp(x,y,z)-vs(x,y,z)*vs(x,y,z)))

        const double ur = (2/15)*3.14159*GravConst(x,y,z)*density(x,y,z)*(fac1/fac1_2)*(fac2*OuterRad(x,y,z)*OuterRad(x,y,z)*r-fac3*r*r*r)

        //put back in cartesian..


        int sgn(double d){ 
            return d<-eps?-1:d>eps;
        }
        
        
        theta = acos(z/(r))
        phi= acos(x/(x*x+y*y))*sgn(y)

        z_disp = cos(theta)*ur
        return z_disp;
    } // disp_z


    // potential
    static double potential(const double x,
                            const double y,
                            const double z) {

        const double r = pow(pow(x,2)+pow(x,2)+pow(x,2),0.5);

        return (4/3)*3.14159*GravConst(x,y,z)*density(x,y,z)*r;
    } // potential


    static PetscErrorCode solnkernel_potential(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* s,
                                              void* context) {
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
                                          PetscScalar* s,
                                          void* context) {
        assert(3 == spaceDim);
        assert(3 == numComponents);
        assert(s);

        s[0] = disp_x(x[0], x[1], x[2]);
        s[1] = disp_y(x[0], x[1], x[2]);
        s[2] = disp_z(x[0], x[1], x[2]);

        return 0;
    } // solnkernel_disp

public:

    static
    TestSelfGrav_Data* createData(void) {
        TestSelfGrav_Data* data = new TestSelfGrav_Data();assert(data);

        data->journalName = "SelfGrav3D";

        data->isJacobianLinear = true;

        data->meshFilename = ":UNKNOWN:"; // Set in child class.
        data->boundaryLabel = "boundary";

        data->normalizer.setLengthScale(LENGTHSCALE);
        data->normalizer.setTimeScale(TIMESCALE);
        data->normalizer.setPRESSURESCALE(PRESSURESCALE);
        data->normalizer.computeDensityScale();

        // solnDiscretizations set in derived class.

        // Material information
        data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = { // order must match order of subfields in auxiliary field
            "density",
            "body_force",
            "shear_modulus",
            "bulk_modulus",
        };
        data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // body_force
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

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

        { // potential
            static const PylithInt constrainedDOF[1] = {0};
            static const PylithInt numConstrainedDOF = 1;
            pylith::bc::DirichletUserFn* bc = new pylith::bc::DirichletUserFn();assert(bc);
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

        data->exactSolnFns = const_cast<pylith::testing::MMSTest::solution_fn*>(_exactSolnFns);
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
pylith::TestSelfGrav_Data*
pylith::SelfGrav3D::TetP2(void) {
    TestSelfGrav_Data* data = pylith::_SelfGrav3D::createData();assert(data);

    data->meshFilename = "data/tet.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // disp
        pylith::topology::Field::Discretization(1, 1), // potential
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // body_force
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TetP2


// ------------------------------------------------------------------------------------------------
pylith::TestSelfGrav_Data*
pylith::SelfGrav3D::TetP3(void) {
    TestSelfGrav_Data* data = pylith::_SelfGrav3D::createData();assert(data);

    data->meshFilename = "data/tet.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // disp
        pylith::topology::Field::Discretization(1, 2), // potential
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // body_force
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // TetP3


// ------------------------------------------------------------------------------------------------
pylith::TestSelfGrav_Data*
pylith::SelfGrav3D::HexQ2(void) {
    TestSelfGrav_Data* data = pylith::_SelfGrav3D::createData();assert(data);

    data->meshFilename = "data/hex.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // potential
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // body_force
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // HexQ2


// ------------------------------------------------------------------------------------------------
pylith::TestSelfGrav_Data*
pylith::SelfGrav3D::HexQ3(void) {
    TestSelfGrav_Data* data = pylith::_SelfGrav3D::createData();assert(data);

    data->meshFilename = "data/hex.mesh";

    data->numSolnSubfields = 2;
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(3, 3), // disp
        pylith::topology::Field::Discretization(2, 3), // potential
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
        pylith::topology::Field::Discretization(0, 3), // density
        pylith::topology::Field::Discretization(0, 3), // body_force
        pylith::topology::Field::Discretization(0, 3), // shear_modulus
        pylith::topology::Field::Discretization(0, 3), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    return data;
} // HexQ3


// End of file
