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
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file mmstets/elasticity/TestFaultKin.hh
 *
 * @brief C++ abstract base class for testing faults with prescribed slip.
 */

#if !defined(pylith_mmstests_testfaultkin_hh)
#define pylith_mmstests_testfaultkin_hh

#include "pylith/testing/MMSTest.hh" // ISA MMSTEST

#include "pylith/faults/faultsfwd.hh" // HOLDSA FaultCohesiveKin
#include "pylith/materials/materialsfwd.hh" // HOLDSA Material
#include "pylith/bc/bcfwd.hh" // USES DirichletUserFn
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

namespace pylith {
    namespace mmstests {
        class TestFaultKin;

        class TestFaultKin_Data; // test data
    } // mmstets
} // pylith

class pylith::mmstests::TestFaultKin : public pylith::testing::MMSTest {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    virtual
    void setUp(void);

    /// Deallocate testing data.
    virtual
    void tearDown(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Initialize objects for test.
    void _initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::faults::FaultCohesiveKin* _fault; ///< Fault test subject.
    pylith::materials::Elasticity* _material; ///< Elastic material.
    pylith::bc::DirichletUserFn* _bc; ///< Dirichlet boundary condition.
    TestFaultKin_Data* _data; ///< Test parameters.

}; // class TestFaultKin

// =====================================================================================================================
class pylith::mmstests::TestFaultKin_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestFaultKin_Data(void);

    /// Destructor
    ~TestFaultKin_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    int spaceDim; ///< Spatial dimension of domain.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* boundaryLabel; ///< Group defining domain boundary.

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    PylithReal startTime; ///< Start time for simulation.
    PylithReal endTime; ///< Total time in simulation.
    PylithReal timeStep; ///< Time step in simulation.
    bool isExplicit; ///< True for explicit time stepping.

    int numSolnSubfields; ///< Number of solution fields.
    pylith::topology::Field::Discretization* solnDiscretizations; ///< Discretizations for solution fields.

    pylith::materials::RheologyElasticity* rheology; ///< Elastic rheology for material.
    int matNumAuxSubfields; ///< Number of material auxiliary subfields.
    const char** matAuxSubfields; ///< Names of material auxiliary subfields.
    pylith::topology::Field::Discretization* matAuxDiscretizations; ///< Discretizations for material aux subfields.
    spatialdata::spatialdb::UserFunctionDB* matAuxDB; ///< Spatial database for material auxiliary field.

    pylith::faults::KinSrc* kinsrc; ///< Kinematic description of fault rupture.
    int faultNumAuxSubfields; ///< Number of fault auxiliary subfields.
    const char** faultAuxSubfields; ///< Names of fault auxiliary subfields.
    pylith::topology::Field::Discretization* faultAuxDiscretizations; ///< Discretizations for fault aux subfields.
    spatialdata::spatialdb::UserFunctionDB* faultAuxDB; ///< Spatial database for fault auxiliary field.

}; // TestFaultKin_Data

#endif // pylith_mmstests_testfaultkin_hh

// End of file
