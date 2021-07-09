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
 * @file tests/libtests/feassemble/TestIntegratorDomain.hh
 *
 * @brief C++ class for testing IntegratorDomain.
 */

#if !defined(pylith_feassemble_testintegratordomain_hh)
#define pylith_feassemble_testintegratordomain_hh

#include <cppunit/extensions/HelperMacros.h>
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ResidualKernels
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

/// Namespace for pylith package
namespace pylith {
    namespace feassemble {
        class TestIntegratorDomain;

        class TestIntegratorDomain_Data; // test data
    } // feassemble
} // pylith

class pylith::feassemble::TestIntegratorDomain : public CppUnit::TestFixture, public pylith::utils::GenericComponent {
    // CPPUNIT TEST SUITE //////////////////////////////////////////////////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestIntegratorDomain);

    CPPUNIT_TEST(testAccessors);

    CPPUNIT_TEST(testSetKernels);
    CPPUNIT_TEST(testInitialize);
    CPPUNIT_TEST(testPoststep);
    CPPUNIT_TEST(testUpdateState);
    CPPUNIT_TEST(testComputeResidual);
    CPPUNIT_TEST(testComputeJacobian);

    CPPUNIT_TEST_SUITE_END_ABSTRACT();

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    virtual
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test getPhysicsDomainMesh(), getAuxiliaryField(), getDerivedField(), getMaterialId(), setMaterialId().
    void testAccessors(void);

    /// Test setKernelsRHSResidual(), setKernelsLHSResidual(), setKernelsLHSJacobian().
    void testSetKernels(void);

    /// Test initialize().
    void testInitialize(void);

    /// Test poststep().
    void testPoststep(void);

    /// Test updateState().
    void testUpdateState(void);

    /// Test computeRHSResidual(), computeLHSResidual().
    void testComputeResidual(void);

    /// Test computeLHSJacobian().
    void testComputeJacobian(void);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Do minimal initilaization of test data.
    void _initializeMin(void);

    /// Do full initilaization of test data.
    void _initializeFull(void);

    /// Setup and populate solution fields.
    virtual
    void _setupSolutionFields(void) = 0;

    /** Set field to zero on the boundary.
     *
     * @param[out] field Field in which to set boundary values to zero.
     */
    void _zeroBoundary(pylith::topology::Field* field);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::feassemble::IntegratorDomain* _integrator; ///< Test subject.
    TestIntegratorDomain_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Fields* _solutionFields; ///< Container for solution fields.

}; // class TestIntegratorDomain

// =====================================================================================================================
class pylith::feassemble::TestIntegratorDomain_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestIntegratorDomain_Data(void);

    /// Destructor
    ~TestIntegratorDomain_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    int dimension; ///< Dimension of domain integrator.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    const char* boundaryLabel; ///< Label for group defining domain boundary.
    PylithInt materialId; ///< Identifier of cells in integration domain.

    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    PylithReal t; ///< Time for solution in simulation.
    PylithReal dt; ///< Time step in simulation.
    PylithInt tindex; ///< Time step index in simulation.
    PylithReal s_tshift; ///< Time shift for LHS Jacobian.

    int numSolutionSubfields; ///< Number of solution fields.
    pylith::topology::Field::Discretization* solutionDiscretizations; ///< Discretizations for solution fields.
    spatialdata::spatialdb::UserFunctionDB* solutionDB; ///< Spatial database with solution.
    spatialdata::spatialdb::UserFunctionDB* perturbationDB; ///< Spatial database with solution + perturbation.

    int numAuxiliarySubfields; ///< Number of auxiliary subfields.
    const char** auxiliarySubfields; ///< Names of auxiliary subfields.
    pylith::topology::Field::Discretization* auxiliaryDiscretizations; ///< Discretizations for auxiliary subfields.
    spatialdata::spatialdb::UserFunctionDB* auxiliaryDB; ///< Spatial database with auxiliary field.
    spatialdata::spatialdb::UserFunctionDB* auxiliaryUpdateDB; ///< Spatial database with updated auxiliary field.

    std::vector<IntegratorDomain::ResidualKernels> kernelsRHSResidual;
    std::vector<IntegratorDomain::ResidualKernels> kernelsLHSResidual;
    std::vector<IntegratorDomain::JacobianKernels> kernelsLHSJacobian;
    std::vector<IntegratorDomain::ProjectKernels> kernelsUpdateStateVars;
    std::vector<IntegratorDomain::ProjectKernels> kernelsDerivedField;
    bool hasLHSJacobianLumpedInv;
};

#endif // pylith_feassemble_testintegratordomain_hh

// End of file
