// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/faults/TestFaultCohesiveDynL.hh
 *
 * @brief C++ TestFaultCohesiveDynL object
 *
 * C++ unit testing for FaultCohesiveDynL.
 */

#if !defined(pylith_faults_testfaultcohesivedynl_hh)
#define pylith_faults_testfaultcohesivedynl_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/faults/faultsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // USES Mesh, SubMesh
#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Quadrature
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include <vector> // HASA std::vector
/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestFaultCohesiveDynL;

    class CohesiveDynLData;
  } // faults
} // pylith

/// C++ unit testing for FaultCohesiveDynL
class pylith::faults::TestFaultCohesiveDynL: public CppUnit::TestFixture { // class TestFaultCohesiveDynL

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
CPPUNIT_TEST_SUITE( TestFaultCohesiveDynL );

    CPPUNIT_TEST( testConstructor );
    CPPUNIT_TEST( testDBInitialTract );

    // Tests in derived classes:
    // testInitialize()
    // testConstrainSolnSpaceStick()
    // testConstrainSolnSpaceSlip()
    // testConstrainSolnSpaceOpen()
    // testUpdateStateVars()
    // testCalcTractions()

  CPPUNIT_TEST_SUITE_END();

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

  CohesiveDynLData* _data; ///< Data for testing
  feassemble::Quadrature<topology::SubMesh>* _quadrature; ///< Fault quad.
  spatialdata::spatialdb::SpatialDB* _dbInitialTract; ///< Initial tractions.
  bool _flipFault; ///< If true, flip fault orientation.

  // PUBLIC METHODS /////////////////////////////////////////////////////
public:

  /// Setup testing data.
  void setUp(void);

  /// Tear down testing data.
  void tearDown(void);

  /// Test constructor.
  void testConstructor(void);

  /// Test dbInitialTract().
  void testDBInitialTract(void);

  /// Test initialize().
  void testInitialize(void);

  /// Test constrainSolnSpace() for sticking case.
  void testConstrainSolnSpaceStick(void);

  /// Test constrainSolnSpace() for slipping case.
  void testConstrainSolnSpaceSlip(void);

  /// Test constrainSolnSpace for fault opening case().
  void testConstrainSolnSpaceOpen(void);

  /// Test updateStateVars().
  void testUpdateStateVars(void);

  /// Test _calcTractions().
  void testCalcTractions(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private:

  /** Initialize FaultCohesiveDynL interface condition.
   *
   * @param mesh PETSc mesh to initialize
   * @param fault Cohesive fault interface condition to initialize.
   * @param fields Solution fields.
   */
  void _initialize(topology::Mesh* const mesh,
      FaultCohesiveDynL* const fault,
      topology::SolutionFields* const fields);

  /** Set values for fields and Jacobian.
   *
   * @pre Must call _initialize() before _setFieldsJacobian to set solution
   * field. Note: Call to Jacobian constructor should be after call to
   * _initialize() so that the solution field is setup.
   *
   * @param mesh PETSc mesh to initialize
   * @param fault Cohesive fault interface condition to initialize.
   * @param fields Solution fields.
   * @param jacobian Jacobian sparse matrix.
   * @param fieldIncrVals Values for solution increment field.
   */
  void _setFieldsJacobian(topology::Mesh* const mesh,
      FaultCohesiveDynL* const fault,
      topology::SolutionFields* const fields,
      topology::Jacobian* const jacobian,
      const double* const fieldIncrVals);

  /** Determine if vertex is a Lagrange multiplier constraint vertex.
   *
   * @param vertex Label of vertex.
   *
   * @returns True if vertex is a constraint vertex, false otherwise.
   */
  bool _isConstraintVertex(const int vertex) const;

}; // class TestFaultCohesiveDynL

#endif // pylith_faults_testfaultcohesivedynl_hh

// End of file 
