# C++ unit tests

The C++ unit tests target verification at the scale of individual C++ class methods.
This isolates bugs very close to their origin, so they are highly effective at catching simple errors, especially memory errors such as invalid reads, writes, and leaks.
We use the [CppUnit](https://www.freedesktop.org/wiki/Software/cppunit/) testing framework for constructing and running tests.

For a given C++ class in `libsrc/pylith`, we create a separate test class in `tests/libtests`.
For example, we create `pylith::problems::TestPhysics` in `tests/libtsts/problems/TestPhysics.*` to test the `pylith::problems::Physics` class in `libsrc/pylith/problems/Physics.*`.
For simple classes that can be fully tested with a single test case, we put the class declaration and implementation in a single `.cc` file.
For classes that require testing with multiple alternative test cases, we usually define both the test class and a data class in a header file (for example `TestSolutionFactory.hh`).
The classes are implemented in a corresponding `.cc` file (for example, `TestSolutionFactory.cc`) with the test cases defined and implemented in a `_Cases.cc` file (for example, `TestSolutionFactory_Cases.cc`).

The test class inherits from `CppUnit::TestFixture` and defines a list of test methods.
The class methods `setUp()` and `tearDown()` are run before and after each test, respectively.
If a class requires significant initialization, we usually put that in protected methods.

## CppUnit macros

* Generating the test suite
  * **`CPPUNIT_TEST_SUITE(ClassName);`** Create a test suite for class `ClassName`.
  * **`CPPUNIT_TEST(classMethod);`** Add `classMethod` as a test in the test suite.
  * **`CPPUNIT_TEST_SUITE_END();`** Mark end of test suite.
* Implementing tests
  * **`CPPUNIT_ASSERT(condition)`** Test that `condition==true`.
  * **`CPPUNIT_ASSERT_EQUAL(valueExpected, valueTest);`** Test that `valueExpected == valueTest`.
  * **`CPPUNIT_ASSERT_DOUBLES_EQUAL(valueExpected, valueTest, tolerance);`** Test that `valueExpected == valueTest within tolerance`.
  * **`CPPUNIT_FAIL(msg);`** Force test failure and include msg string in test failure message.
  * **`CPPUNIT_ASSERT_MESSAGE(msg, condition);`** Include msg string in test failure message.
  * **`CPPUNIT_ASSERT_EQUAL_MESSAGE(msg, valueExpected, valueTest);`**
  * **`CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(msg, valueExpected, valueTest, tolernace);`**

Within the test implementations, we use `CPPUNIT_ASSERT(condition)` where we would normally use `assert(condition)`.
This results in the corresponding test failing and the remaining tests are executed.
If an `assert(condition)` fails, the test driver aborts and the remaining tests are not executed.

```{code-block} c++
---
caption: Illustration of a C++ header file for a test class without a separate data class.
---
#if !defined(pylith_problems_testphysics_hh)
#define pylith_problems_testphysics_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/problems/problemsfwd.hh" // HOLDSA Physics
#include "pylith/topology/topologyfwd.hh" // HOLDSA Mesh, Field

// Declaration of class within the pylith namespace.
namespace pylith {
    namespace problems {
        class TestPhysics;
    } // problems
} // pylith

class pylith::problems::TestPhysics : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(TestPhysics); // CppUnit macro used to define the `TestPhysics` test suite.

    CPPUNIT_TEST(testSetNormalizer); // CppUnit macro to add test implemented by class method.
    CPPUNIT_TEST(testSetAuxiliaryFieldDB);
    CPPUNIT_TEST(testSetAuxiliarySubfieldDiscretization);
    CPPUNIT_TEST(testObservers);
    CPPUNIT_TEST(testGetKernelConstants);
    CPPUNIT_TEST(testVerifyConfiguration);
    CPPUNIT_TEST(testCreateIntegrator);
    CPPUNIT_TEST(testCreateConstraint);
    CPPUNIT_TEST(testCreateAuxiliaryField);
    CPPUNIT_TEST(testCreateDerivedField);

    CPPUNIT_TEST_SUITE_END(); // End of test suite declaration.

public:

    // Methods to setup and clean up tests.
    void setUp(void);
    void tearDown(void);

    // Methods that implement tests. The naming convention is testMethodName.
    // All test methods must return void and not have any arguments.
    void testSetNormalizer(void);
    void testSetAuxiliaryFieldDB(void);
    void testSetAuxiliarySubfieldDiscretization(void);
    void testObservers(void);
    void testGetKernelConstants(void);
    void testVerifyConfiguration(void);
    void testCreateIntegrator(void);
    void testCreateConstraint(void);
    void testCreateAuxiliaryField(void);
    void testCreateDerivedField(void);

private:

    pylith::problems::Physics* _physics; // Test subject.
    pylith::topology::Mesh* _mesh; // Mesh for test subject.
    pylith::topology::Field* _solution; // Solution field for test subject.

}; // class TestPhysics

#endif // pylith_problems_testphysics_hh
```

```{code-block} c++
---
caption: Illustration of a C++ implementation file for a test class.
---
// Instantiate the test suite.
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::problems::TestPhysics);

void
pylith::problems::TestPhysics::setUp(void) {
    // We use CPPUNIT_ASSERT() where we normally would use assert().
    _physics = new PhysicsStub();CPPUNIT_ASSERT(_physics);

    _mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    _solution = new pylith::topology::Field(*_mesh);CPPUNIT_ASSERT(_solution);
} // setUp


void
pylith::problems::TestPhysics::tearDown(void) {
    delete _physics;_physics = NULL;
    delete _solution;_solution = NULL;
    delete _mesh;_mesh = NULL;
} // tearDown


void
pylith::problems::TestPhysics::testSetNormalizer(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const PylithReal lengthScale = 3.0;
    normalizer.setLengthScale(lengthScale);

    CPPUNIT_ASSERT(_physics);
    _physics->setNormalizer(normalizer);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(lengthScale, _physics->_normalizer->getLengthScale(), 1.0e-6);

    PYLITH_METHOD_END;
} // testSetNormalizer


void
pylith::problems::TestPhysics::testSetAuxiliaryFieldDB(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::spatialdb::UniformDB db;
    db.setLabel("test db");

    CPPUNIT_ASSERT(_physics);
    _physics->setAuxiliaryFieldDB(&db);

    const pylith::feassemble::AuxiliaryFactory* factory = _physics->_getAuxiliaryFactory();CPPUNIT_ASSERT(factory);
    const spatialdata::spatialdb::SpatialDB* queryDB = factory->getQueryDB();CPPUNIT_ASSERT(queryDB);
    CPPUNIT_ASSERT_EQUAL(db.getLabel(),queryDB->getLabel());

    PYLITH_METHOD_END;
} // testSetAuxiliaryFieldDB

void
pylith::problems::TestPhysics::testGetKernelConstants(void) {
    PYLITH_METHOD_BEGIN;

    const size_t numConstants = 2;
    const PylithReal constantsE[numConstants] = { -1.1, 4.4 };

    CPPUNIT_ASSERT(_physics);
    _physics->_kernelConstants = pylith::real_array(constantsE, numConstants);

    const PylithReal dt = 2.0;
    const pylith::real_array& constants = _physics->getKernelConstants(dt);

    CPPUNIT_ASSERT_EQUAL(numConstants, constants.size());
    for (size_t i = 0; i < numConstants; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(constantsE[i], constants[i], 1.0e-6);
    } // for

    PYLITH_METHOD_END;
} // testGetKernelConstants

```
