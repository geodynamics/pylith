# MMS Test Redesign: Approach C — Detailed Design Specification

## 1. Overview

This document specifies the redesign of MMS (Method of Manufactured Solutions) test
infrastructure in PyLith. The design uses an **Equation Mixin + Shared TestCase Base**
pattern to eliminate code duplication across governing equations while keeping each
test case simple and self-contained.

### 1.1 Goals

- Make **minimal changes** to the existing `MMSTest` base class (`tests/src/MMSTest.hh`).
- Eliminate the per-equation `Test*` classes and `*_Data` structs that are currently
  duplicated across directories.
- Keep individual test case classes simple: define analytical functions, configure
  parameters, expose mesh/discretization variants.
- Support adding new governing equations (thermoplasticity, poroelasticity with
  viscoelastic rheologies) without modifying existing shared infrastructure.
- Support composing faults orthogonally with any equation type.
- All auxiliary fields use `spatialdata::spatialdb::UserFunctionDB` (analytical functions).
- All solution fields use PETSc kernel functions (`solution_fn` signature).

### 1.2 Constraints

- `MMSTest` (`tests/src/MMSTest.hh` and `tests/src/MMSTest.cc`) remains unchanged.
- Body force and gravity continue as special-case fields on the problem.
- Neumann BCs remain a separate BC class (`NeumannUserFn`).
- Single kinematic source per fault for now; design must not preclude multiple sources.
- 3D fault directories are planned.
- Additional governing equations are planned: thermoplasticity, poroelasticity with
  viscoelastic bulk rheologies (5-7+ equation/rheology combinations total).

### 1.3 Assumptions

1. The coordinate system is always Cartesian (`spatialdata::geocoords::CSCart`).
2. Scales follow elasticity conventions via `pylith::scales::ElasticityScales`.
3. Each `tests/mmstests/<equation>/<variant>` directory uses a fixed set of mesh files
   in its `data/` subdirectory. Different test cases within a directory share the mesh
   but vary the analytical solution.
4. Multiple cell types (Tri, Quad, Tet, Hex) and polynomial orders (P1-P4, Q1-Q4) are
   tested for each analytical solution.
5. Some test cases omit certain MMS checks (e.g., skip `testJacobianFiniteDiff` or
   Jacobian tests entirely for dynamic cases).
6. The test framework is Catch2 with the existing `driver_catch2.cc` entry point.
7. Build system uses GNU Autotools (`Makefile.am`).

---

## 2. Class Hierarchy

```
pylith::testing::MMSTest           (UNCHANGED — tests/src/MMSTest.hh)
    │
    └── pylith::testing::MMSTestCase   (NEW — tests/src/MMSTestCase.hh)
            │
            │  Equation Mixins (one per governing equation):
            ├── pylith::testing::ElasticitySetup
            ├── pylith::testing::IncompressibleElasticitySetup
            ├── pylith::testing::PoroelasticitySetup
            ├── pylith::testing::ThermoplasticitySetup          (future)
            │
            │  Fault Mixin (orthogonal to equation type):
            └── pylith::testing::FaultKinSetup
```

Concrete test cases compose via multiple inheritance:

```cpp
// No-fault linear elasticity 2D
class UniformStrain2D : public MMSTestCase, public ElasticitySetup { ... };

// Fault + linear elasticity 2D
class TwoBlocksStatic : public MMSTestCase, public ElasticitySetup, public FaultKinSetup { ... };

// Poroelasticity (no fault)
class PressureGradient : public MMSTestCase, public PoroelasticitySetup { ... };

// Incompressible elasticity
class UniformShear2D : public MMSTestCase, public IncompressibleElasticitySetup { ... };

// Future: poroelasticity + fault 3D
class PoroFault3D : public MMSTestCase, public PoroelasticitySetup, public FaultKinSetup { ... };
```

---

## 3. `MMSTestCase` — The Shared Base Class

### 3.1 Responsibilities

- Owns all data common to every MMS test regardless of equation type.
- Implements the fixed initialization hook sequence.
- Provides default implementations for mesh I/O, nondimensionalization, problem
  wiring, and exact solution registration.
- Declares pure virtual hooks that equation mixins must implement.
- Declares virtual hooks with no-op defaults that the fault mixin overrides.

### 3.2 Header: `tests/src/MMSTestCase.hh`

```cpp
#pragma once

#include "tests/src/MMSTest.hh"

#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "spatialdata/geocoords/CSCart.hh"
#include "spatialdata/spatialdb/spatialdbfwd.hh"
#include "pylith/scales/Scales.hh"
#include "pylith/problems/Physics.hh"
#include "pylith/topology/Field.hh"
#include "pylith/bc/bcfwd.hh"
#include "pylith/problems/problemsfwd.hh"

#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class MMSTestCase;
    }
}

class pylith::testing::MMSTestCase : public pylith::testing::MMSTest {
    // PUBLIC TYPES ///////////////////////////////////////////////////////////////////////////
public:

    /// Bitmask for which MMS checks to run.
    enum TestCheck {
        CHECK_DISCRETIZATION  = 0x1,
        CHECK_RESIDUAL        = 0x2,
        CHECK_JACOBIAN_TAYLOR = 0x4,
        CHECK_JACOBIAN_FD     = 0x8,
        CHECK_ALL = CHECK_DISCRETIZATION | CHECK_RESIDUAL
                    | CHECK_JACOBIAN_TAYLOR | CHECK_JACOBIAN_FD,
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////
public:

    MMSTestCase(void);
    virtual ~MMSTestCase(void);

    /** Query whether a specific check is enabled.
     * @param[in] check Bitmask value from TestCheck enum.
     * @returns true if the check is enabled.
     */
    bool checkEnabled(unsigned int check) const;

    // PROTECTED METHODS — HOOK SEQUENCE //////////////////////////////////////////////////////
protected:

    /** Full initialization. Calls hooks in fixed order.
     * Overrides MMSTest::_initialize().
     */
    void _initialize(void) override;

    /** Register exact solution with PetscDS.
     * Overrides MMSTest::_setExactSolution().
     */
    void _setExactSolution(void) override;

    // --- Hooks (virtual, override in mixins or concrete classes) ---

    /// Read mesh from file. Default handles MeshIOAscii and MeshIOPetsc.
    virtual void _readMesh(void);

    /// Transform mesh topology for faults. Default: no-op.
    virtual void _setupFaultTopology(void);

    /// Nondimensionalize the mesh. Default uses _cs and _scales.
    virtual void _nondimensionalize(void);

    /// Create and configure material objects. PURE VIRTUAL — equation mixin provides.
    virtual void _setupMaterials(void) = 0;

    /** Register solution subfields with SolutionFactory.
     * PURE VIRTUAL — equation mixin provides.
     * @param[in] factory SolutionFactory for the solution field.
     */
    virtual void _setupSolutionSubfields(
        pylith::problems::SolutionFactory& factory) = 0;

    /// Configure fault auxiliary fields and kinematic sources. Default: no-op.
    virtual void _setupFaultFields(void);

    /// Assemble the problem: set materials, BCs, faults, time config. Default provided.
    virtual void _setupProblem(void);

    /// Register exact solution on domain PetscDS. Default iterates _exactSolnFns.
    virtual void _registerExactSolution(PetscDS ds, PetscDM dm);

    /// Register exact solution on fault cohesive cells. Default: no-op.
    virtual void _registerExactSolutionFault(PetscDM dm);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////
protected:

    // --- Mesh ---
    std::string _meshFilename;
    std::string _meshOptions;
    bool _useAsciiMesh;

    // --- Coordinate system and scales ---
    spatialdata::geocoords::CSCart _cs;
    pylith::scales::Scales _scales;
    int _spaceDim;

    // --- Time stepping ---
    PylithReal _t;
    PylithReal _dt;
    pylith::problems::Physics::FormulationEnum _formulation;

    // --- Gravity (optional, nullptr if not used) ---
    spatialdata::spatialdb::GravityField* _gravityField;

    // --- Boundary conditions (owned) ---
    std::vector<pylith::bc::BoundaryCondition*> _bcs;

    // --- Solution discretizations and exact solution functions ---
    std::vector<pylith::topology::Field::Discretization> _solnDiscretizations;
    std::vector<solution_fn> _exactSolnFns;
    std::vector<solution_fn> _exactSolnDotFns;

    // --- Test configuration ---
    unsigned int _enabledChecks;  ///< Bitmask of TestCheck values.

}; // class MMSTestCase
```

### 3.3 Implementation: `tests/src/MMSTestCase.cc`

#### Constructor

```cpp
pylith::testing::MMSTestCase::MMSTestCase(void) :
    _useAsciiMesh(true),
    _spaceDim(2),
    _t(0.0),
    _dt(0.05),
    _formulation(pylith::problems::Physics::QUASISTATIC),
    _gravityField(nullptr),
    _enabledChecks(CHECK_ALL) {}
```

#### Destructor

```cpp
pylith::testing::MMSTestCase::~MMSTestCase(void) {
    for (auto* bc : _bcs) { delete bc; }
    _bcs.clear();
    delete _gravityField; _gravityField = nullptr;
}
```

#### `checkEnabled`

```cpp
bool
pylith::testing::MMSTestCase::checkEnabled(unsigned int check) const {
    return (_enabledChecks & check) != 0;
}
```

#### `_initialize` — The Hook Sequence

```cpp
void
pylith::testing::MMSTestCase::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(_mesh);

    // 1. Read mesh
    _readMesh();

    // 2. Insert cohesive cells for faults (no-op if no faults)
    _setupFaultTopology();

    // 3. Nondimensionalize
    _nondimensionalize();

    // 4. Create and configure materials
    _setupMaterials();

    // 5. Configure fault auxiliary fields (no-op if no faults)
    _setupFaultFields();

    // 6. Wire up the problem object
    _setupProblem();

    // 7. Create solution field and register subfields
    REQUIRE(!_solution);
    _solution = new pylith::topology::Field(*_mesh);
    REQUIRE(_solution);
    _solution->setLabel("solution");
    pylith::problems::SolutionFactory factory(*_solution, _scales);
    _setupSolutionSubfields(factory);
    _problem->setSolution(_solution);

    // 8. Delegate to MMSTest base for PETSc TS/SNES setup
    pylith::testing::MMSTest::_initialize();

    PYLITH_METHOD_END;
}
```

#### `_readMesh`

```cpp
void
pylith::testing::MMSTestCase::_readMesh(void) {
    REQUIRE(!_meshFilename.empty());

    if (_useAsciiMesh) {
        pylith::meshio::MeshIOAscii iohandler;
        iohandler.setFilename(_meshFilename.c_str());
        iohandler.read(_mesh);
    } else {
        if (!_meshOptions.empty()) {
            PylithCallPetsc(PetscOptionsInsertString(nullptr, _meshOptions.c_str()));
        }
        pylith::meshio::MeshIOPetsc iohandler;
        iohandler.setFilename(_meshFilename.c_str());
        iohandler.read(_mesh);
    }
    REQUIRE(_mesh);
    REQUIRE(pylith::topology::MeshOps::getNumCells(*_mesh) > 0);
    REQUIRE(pylith::topology::MeshOps::getNumVertices(*_mesh) > 0);
}
```

#### `_setupFaultTopology` (default no-op)

```cpp
void pylith::testing::MMSTestCase::_setupFaultTopology(void) {}
```

#### `_nondimensionalize`

```cpp
void
pylith::testing::MMSTestCase::_nondimensionalize(void) {
    _mesh->setCoordSys(&_cs);
    pylith::topology::MeshOps::nondimensionalize(_mesh, _scales);
}
```

#### `_setupFaultFields` (default no-op)

```cpp
void pylith::testing::MMSTestCase::_setupFaultFields(void) {}
```

#### `_setupProblem`

This is the default implementation. Subclasses may override if they need
custom problem configuration. The `_getMaterials()`, `_getFaults()` helper
methods are provided by the mixins (see Section 4 and 6).

```cpp
void
pylith::testing::MMSTestCase::_setupProblem(void) {
    REQUIRE(_problem);

    _problem->setScales(_scales);
    _problem->setGravityField(_gravityField);

    // Materials — provided by equation mixin via virtual _getMaterials()
    auto materials = _getMaterialPtrs();
    _problem->setMaterials(materials.data(), materials.size());

    // Faults — provided by fault mixin via virtual _getFaultPtrs()
    auto faults = _getFaultPtrs();
    if (!faults.empty()) {
        _problem->setInterfaces(faults.data(), faults.size());
    }

    // Boundary conditions
    _problem->setBoundaryConditions(_bcs.data(), _bcs.size());

    // Time configuration
    _problem->setStartTime(_t);
    _problem->setEndTime(_t + _dt);
    _problem->setInitialTimeStep(_dt);
    _problem->setFormulation(_formulation);
}
```

**Note:** `_getMaterialPtrs()` and `_getFaultPtrs()` are virtual methods declared
in `MMSTestCase` with default implementations returning empty vectors. Equation
mixins and `FaultKinSetup` override them. See Sections 4 and 6.

```cpp
// In MMSTestCase.hh, protected section:
virtual std::vector<pylith::materials::Material*> _getMaterialPtrs(void);
virtual std::vector<pylith::faults::FaultCohesive*> _getFaultPtrs(void);
```

Default implementations:
```cpp
std::vector<pylith::materials::Material*>
pylith::testing::MMSTestCase::_getMaterialPtrs(void) { return {}; }

std::vector<pylith::faults::FaultCohesive*>
pylith::testing::MMSTestCase::_getFaultPtrs(void) { return {}; }
```

#### `_setExactSolution`

```cpp
void
pylith::testing::MMSTestCase::_setExactSolution(void) {
    REQUIRE(!_exactSolnFns.empty());

    const pylith::topology::Field* solution = _problem->getSolution();
    REQUIRE(solution);

    PetscDM dm = solution->getDM();
    PetscDS ds = nullptr;
    PylithCallPetsc(DMGetDS(dm, &ds));

    _registerExactSolution(ds, dm);
    _registerExactSolutionFault(dm);
}
```

#### `_registerExactSolution` (default)

```cpp
void
pylith::testing::MMSTestCase::_registerExactSolution(PetscDS ds, PetscDM dm) {
    for (size_t i = 0; i < _exactSolnFns.size(); ++i) {
        PylithCallPetsc(PetscDSSetExactSolution(ds, i, _exactSolnFns[i], nullptr));
        if (i < _exactSolnDotFns.size() && _exactSolnDotFns[i]) {
            PylithCallPetsc(PetscDSSetExactSolutionTimeDerivative(
                ds, i, _exactSolnDotFns[i], nullptr));
        }
    }
}
```

#### `_registerExactSolutionFault` (default no-op)

```cpp
void pylith::testing::MMSTestCase::_registerExactSolutionFault(PetscDM dm) {}
```

---

## 4. Equation Mixins

Each equation mixin is a non-polymorphic helper class that:
- Owns the material object(s) and rheology for its equation type.
- Owns the `UserFunctionDB` for material auxiliary fields.
- Provides protected helper methods that the concrete test case calls from its
  hook overrides.
- Provides `_getMaterialPtrs()` override via a helper that the concrete class delegates to.

### 4.1 `ElasticitySetup`

**File:** `tests/src/ElasticitySetup.hh`

```cpp
#pragma once

#include "pylith/materials/Elasticity.hh"
#include "pylith/materials/IsotropicLinearElasticity.hh"
#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "pylith/topology/Field.hh"
#include "pylith/problems/problemsfwd.hh"

#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class ElasticitySetup;
    }
}

class pylith::testing::ElasticitySetup {
protected:

    ElasticitySetup(void);
    virtual ~ElasticitySetup(void);

    /// Configuration for a single material block.
    struct MaterialConfig {
        std::string name;       ///< e.g., "material-id=24"
        int labelValue;         ///< Label value in mesh.
        bool useBodyForce;      ///< Enable body force.
        bool useReferenceState; ///< Enable reference state in rheology.
    };

    // --- Setup helpers (called from hook overrides) ---

    /** Create Elasticity material objects from _matConfigs.
     * Sets auxiliary field DB and discretizations on each material.
     * @param[in] spaceDim Spatial dimension.
     * @param[in] formulation Quasistatic or dynamic.
     */
    void _createElasticityMaterials(int spaceDim,
                                    pylith::problems::Physics::FormulationEnum formulation);

    /** Add displacement (and optionally velocity) to solution.
     * @param[in] factory SolutionFactory reference.
     * @param[in] formulation Quasistatic or dynamic.
     */
    void _addElasticitySolutionSubfields(
        pylith::problems::SolutionFactory& factory,
        pylith::problems::Physics::FormulationEnum formulation);

    /** Get raw material pointers for problem setup.
     * @returns Vector of Material* pointing to owned Elasticity objects.
     */
    std::vector<pylith::materials::Material*> _getElasticityMaterialPtrs(void);

    // --- Data members ---

    std::vector<MaterialConfig> _matConfigs;

    /// Auxiliary field spatial database (analytical functions for properties).
    spatialdata::spatialdb::UserFunctionDB _matAuxDB;

    /// Names of material auxiliary subfields.
    std::vector<std::string> _matAuxSubfields;

    /// Discretizations for material auxiliary subfields.
    std::vector<pylith::topology::Field::Discretization> _matAuxDiscretizations;

    /// Rheology object (shared across all materials in a test).
    pylith::materials::IsotropicLinearElasticity _rheology;

    /// Owned material objects.
    std::vector<pylith::materials::Elasticity*> _elasticityMaterials;

}; // class ElasticitySetup
```

**File:** `tests/src/ElasticitySetup.cc`

```cpp
pylith::testing::ElasticitySetup::ElasticitySetup(void) {
    _matAuxDB.setDescription("material auxiliary field spatial database");
}

pylith::testing::ElasticitySetup::~ElasticitySetup(void) {
    for (auto* mat : _elasticityMaterials) { delete mat; }
    _elasticityMaterials.clear();
}

void
pylith::testing::ElasticitySetup::_createElasticityMaterials(
    int spaceDim,
    pylith::problems::Physics::FormulationEnum formulation) {

    _elasticityMaterials.resize(_matConfigs.size());
    for (size_t i = 0; i < _matConfigs.size(); ++i) {
        const auto& cfg = _matConfigs[i];
        auto* material = new pylith::materials::Elasticity();
        material->setFormulation(formulation);
        material->useBodyForce(cfg.useBodyForce);
        material->setIdentifier("elasticity");
        material->setName(cfg.name.c_str());
        material->setLabelValue(cfg.labelValue);
        material->setBulkRheology(&_rheology);
        _rheology.useReferenceState(cfg.useReferenceState);

        material->setAuxiliaryFieldDB(&_matAuxDB);
        for (size_t j = 0; j < _matAuxSubfields.size(); ++j) {
            const auto& disc = _matAuxDiscretizations[j];
            material->setAuxiliarySubfieldDiscretization(
                _matAuxSubfields[j].c_str(),
                disc.basisOrder, disc.quadOrder, spaceDim,
                pylith::topology::FieldBase::DEFAULT_BASIS,
                disc.feSpace, disc.isBasisContinuous);
        }
        _elasticityMaterials[i] = material;
    }
}

void
pylith::testing::ElasticitySetup::_addElasticitySolutionSubfields(
    pylith::problems::SolutionFactory& factory,
    pylith::problems::Physics::FormulationEnum formulation) {

    // _solnDiscretizations is on MMSTestCase; accessed via concrete class.
    // This method is called by the concrete class which passes the factory.
    // The concrete class accesses its own _solnDiscretizations.
    // Implementation note: the concrete class calls:
    //   factory.addDisplacement(_solnDiscretizations[0]);
    //   if (DYNAMIC) factory.addVelocity(_solnDiscretizations[1]);
    // This helper is a convenience; see Section 7 for concrete usage.
}

std::vector<pylith::materials::Material*>
pylith::testing::ElasticitySetup::_getElasticityMaterialPtrs(void) {
    std::vector<pylith::materials::Material*> ptrs;
    ptrs.reserve(_elasticityMaterials.size());
    for (auto* mat : _elasticityMaterials) {
        ptrs.push_back(mat);
    }
    return ptrs;
}
```

### 4.2 `IncompressibleElasticitySetup`

Same pattern as `ElasticitySetup` but:
- Uses `pylith::materials::IncompressibleElasticity` material type.
- Uses `pylith::materials::IsotropicLinearIncompElasticity` rheology.
- `_addSolutionSubfields` calls `factory.addDisplacement(...)` and
  `factory.addPressure(...)`.
- For dynamic formulation, also adds velocity.

**File:** `tests/src/IncompressibleElasticitySetup.hh`

```cpp
#pragma once

#include "pylith/materials/IncompressibleElasticity.hh"
#include "pylith/materials/IsotropicLinearIncompElasticity.hh"
#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "pylith/topology/Field.hh"

#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class IncompressibleElasticitySetup;
    }
}

class pylith::testing::IncompressibleElasticitySetup {
protected:
    IncompressibleElasticitySetup(void);
    virtual ~IncompressibleElasticitySetup(void);

    void _createIncompressibleMaterials(
        int spaceDim,
        pylith::problems::Physics::FormulationEnum formulation);

    std::vector<pylith::materials::Material*> _getIncompressibleMaterialPtrs(void);

    // Data members — same structure as ElasticitySetup
    struct MaterialConfig {
        std::string name;
        int labelValue;
        bool useBodyForce;
    };
    std::vector<MaterialConfig> _incompMatConfigs;
    spatialdata::spatialdb::UserFunctionDB _incompMatAuxDB;
    std::vector<std::string> _incompMatAuxSubfields;
    std::vector<pylith::topology::Field::Discretization> _incompMatAuxDiscretizations;
    pylith::materials::IsotropicLinearIncompElasticity _incompRheology;
    std::vector<pylith::materials::IncompressibleElasticity*> _incompMaterials;
};
```

### 4.3 `PoroelasticitySetup`

- Uses `pylith::materials::Poroelasticity` material type.
- Uses `pylith::materials::IsotropicLinearPoroelasticity` rheology.
- `_addSolutionSubfields` calls `factory.addDisplacement(...)`,
  `factory.addPressure(...)`, `factory.addTraceStrain(...)`.
- For time-dependent cases with 6 solution subfields, also adds
  `factory.addVelocity(...)`, `factory.addPressureDot(...)`,
  `factory.addTraceStrainDot(...)`.

**File:** `tests/src/PoroelasticitySetup.hh`

```cpp
#pragma once

#include "pylith/materials/Poroelasticity.hh"
#include "pylith/materials/IsotropicLinearPoroelasticity.hh"
#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "pylith/topology/Field.hh"

#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class PoroelasticitySetup;
    }
}

class pylith::testing::PoroelasticitySetup {
protected:
    PoroelasticitySetup(void);
    virtual ~PoroelasticitySetup(void);

    void _createPoroelasticityMaterials(
        int spaceDim,
        pylith::problems::Physics::FormulationEnum formulation);

    std::vector<pylith::materials::Material*> _getPoroelasticityMaterialPtrs(void);

    struct MaterialConfig {
        std::string name;
        int labelValue;
        bool useBodyForce;
        bool useTensorPermeability;
    };
    std::vector<MaterialConfig> _poroMatConfigs;
    spatialdata::spatialdb::UserFunctionDB _poroMatAuxDB;
    std::vector<std::string> _poroMatAuxSubfields;
    std::vector<pylith::topology::Field::Discretization> _poroMatAuxDiscretizations;
    pylith::materials::IsotropicLinearPoroelasticity _poroRheology;
    std::vector<pylith::materials::Poroelasticity*> _poroMaterials;
};
```

### 4.4 Future Equation Mixins

Each new equation type (thermoplasticity, viscoelastic poroelasticity, etc.)
follows the same pattern:

1. Create `<Equation>Setup.hh` and `<Equation>Setup.cc` in `tests/src/`.
2. Define the material type, rheology type, and `MaterialConfig` struct.
3. Implement `_create<Equation>Materials()` and `_get<Equation>MaterialPtrs()`.
4. The concrete test case inherits from `MMSTestCase` + the new mixin.

No changes to `MMSTestCase`, `MMSTest`, or any other mixin are required.

---

## 5. Equation Mixin Implementation Details

### 5.1 Material Creation Pattern

All equation mixins follow the same internal pattern in `_create*Materials()`:

```
for each MaterialConfig:
    1. Allocate material of the correct type.
    2. Set formulation (quasistatic/dynamic).
    3. Set identifier, name, label value.
    4. Set bulk rheology.
    5. Configure equation-specific flags (useBodyForce, useReferenceState, etc.).
    6. Set auxiliary field DB to the mixin's UserFunctionDB.
    7. Set auxiliary subfield discretizations from the mixin's vectors.
    8. Store in the mixin's owned materials vector.
```

### 5.2 Solution Subfield Registration Pattern

Each equation mixin provides a helper or the concrete class directly calls
`SolutionFactory` methods in its `_setupSolutionSubfields()` override:

| Equation | Quasistatic subfields | Dynamic subfields |
|----------|----------------------|-------------------|
| Elasticity | displacement | displacement, velocity |
| IncompressibleElasticity | displacement, pressure | displacement, pressure, velocity |
| Poroelasticity | displacement, pressure, trace_strain | displacement, pressure, trace_strain, velocity, pressure_dot, trace_strain_dot |
| Thermoplasticity (future) | displacement, temperature | displacement, temperature, velocity |

For fault tests, `addLagrangeMultiplierFault(...)` is appended after the domain
subfields (see Section 6).

---

## 6. `FaultKinSetup` Mixin

### 6.1 Responsibilities

- Owns `FaultCohesiveKin` objects and `KinSrc` objects.
- Owns the `UserFunctionDB` for fault auxiliary fields (slip, initiation time).
- Provides topology transformation (inserting cohesive cells).
- Provides fault field configuration (setting `auxFieldDB` on kinematic sources).
- Provides exact solution registration on cohesive cell `PetscDS`.
- Adds Lagrange multiplier subfield(s) to the solution.

### 6.2 Header: `tests/src/FaultKinSetup.hh`

```cpp
#pragma once

#include "pylith/faults/FaultCohesiveKin.hh"
#include "pylith/faults/KinSrcStep.hh"
#include "pylith/topology/topologyfwd.hh"
#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "pylith/topology/Field.hh"
#include "tests/src/MMSTest.hh"  // for solution_fn typedef

#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class FaultKinSetup;
    }
}

class pylith::testing::FaultKinSetup {
protected:

    FaultKinSetup(void);
    virtual ~FaultKinSetup(void);

    /// Per-fault configuration.
    struct FaultConfig {
        int cohesiveLabelValue;          ///< e.g., 100
        std::string surfaceLabelName;    ///< e.g., "fault_xpos_faces"
        int surfaceLabelValue;           ///< Default 1
        std::string buriedEdgesLabelName; ///< Empty if none
        int buriedEdgesLabelValue;       ///< Default 0
    };

    // --- Setup helpers ---

    /** Insert cohesive cells for all configured faults.
     * Modifies mesh in place (deletes old, replaces with new).
     * @param[inout] mesh Pointer to mesh pointer (will be replaced).
     */
    void _transformFaultTopology(pylith::topology::Mesh*& mesh);

    /** Configure kinematic sources: set auxFieldDB on each KinSrc.
     * @param[in] spaceDim Spatial dimension.
     */
    void _configureFaultFields(int spaceDim);

    /** Add Lagrange multiplier subfield(s) to solution.
     * @param[in] factory SolutionFactory reference.
     */
    void _addLagrangeSubfields(pylith::problems::SolutionFactory& factory);

    /** Get fault pointers for problem setup.
     * @returns Vector of FaultCohesive* pointing to owned objects.
     */
    std::vector<pylith::faults::FaultCohesive*> _getFaultKinPtrs(void);

    /** Register exact solution on fault cohesive cells PetscDS.
     *
     * Gets the cohesive cell DS from DM, registers domain subfield
     * exact solutions (without DM context) and fault subfield exact
     * solutions.
     *
     * @param[in] dm Solution DM.
     * @param[in] allExactSolnFns All exact solution functions (domain + fault).
     * @param[in] allExactSolnDotFns All exact solution dot functions.
     * @param[in] numDomainSubfields Number of domain solution subfields.
     */
    void _registerFaultExactSolution(
        PetscDM dm,
        const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnFns,
        const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnDotFns,
        size_t numDomainSubfields);

    // --- Data members ---

    std::vector<FaultConfig> _faultConfigs;

    /** Kinematic sources — outer vector is per-fault, inner vector is
     * per-source within that fault. Currently inner vector has size 1.
     * Expandable to multiple sources later without interface change.
     */
    std::vector<std::vector<pylith::faults::KinSrc*>> _kinSrcs;

    /// Fault auxiliary field spatial database (slip, initiation_time, etc.).
    spatialdata::spatialdb::UserFunctionDB _faultAuxDB;

    /// Fault auxiliary subfield names.
    std::vector<std::string> _faultAuxSubfields;

    /// Fault auxiliary subfield discretizations.
    std::vector<pylith::topology::Field::Discretization> _faultAuxDiscretizations;

    /// Owned FaultCohesiveKin objects.
    std::vector<pylith::faults::FaultCohesiveKin*> _faultObjects;

    /// Solution discretization(s) for Lagrange multiplier subfield(s).
    std::vector<pylith::topology::Field::Discretization> _lagrangeDiscretizations;

    /// Number of fault solution subfields (for exact solution registration).
    size_t _numFaultSolnSubfields;

}; // class FaultKinSetup
```

### 6.3 Implementation: `tests/src/FaultKinSetup.cc`

```cpp
pylith::testing::FaultKinSetup::FaultKinSetup(void) :
    _numFaultSolnSubfields(0) {
    _faultAuxDB.setDescription("fault auxiliary field spatial database");
}

pylith::testing::FaultKinSetup::~FaultKinSetup(void) {
    for (auto* fault : _faultObjects) { delete fault; }
    _faultObjects.clear();
    for (auto& srcs : _kinSrcs) {
        for (auto* src : srcs) { delete src; }
    }
    _kinSrcs.clear();
}

void
pylith::testing::FaultKinSetup::_transformFaultTopology(
    pylith::topology::Mesh*& mesh) {

    REQUIRE(_faultConfigs.size() == _kinSrcs.size());

    _faultObjects.resize(_faultConfigs.size());
    for (size_t i = 0; i < _faultConfigs.size(); ++i) {
        const auto& cfg = _faultConfigs[i];
        auto* fault = new pylith::faults::FaultCohesiveKin();
        fault->setCohesiveLabelValue(cfg.cohesiveLabelValue);
        fault->setSurfaceLabelName(cfg.surfaceLabelName.c_str());
        if (cfg.surfaceLabelValue != 0) {
            fault->setSurfaceLabelValue(cfg.surfaceLabelValue);
        }
        if (!cfg.buriedEdgesLabelName.empty()) {
            fault->setBuriedEdgesLabelName(cfg.buriedEdgesLabelName.c_str());
            fault->setBuriedEdgesLabelValue(cfg.buriedEdgesLabelValue);
        }

        // Set kinematic ruptures
        const auto& srcs = _kinSrcs[i];
        std::vector<const char*> names;
        std::vector<pylith::faults::KinSrc*> srcPtrs;
        for (size_t j = 0; j < srcs.size(); ++j) {
            std::string name = "rupture_" + std::to_string(j);
            // Store name persistently (in a member if needed)
            names.push_back("rupture");  // simplified for single source
            srcPtrs.push_back(srcs[j]);
        }
        fault->setEqRuptures(names.data(), names.size(),
                             srcPtrs.data(), srcPtrs.size());

        // Transform topology
        pylith::topology::Mesh* meshNew = fault->transformTopology(mesh);
        delete mesh;
        mesh = meshNew;

        _faultObjects[i] = fault;
    }
}

void
pylith::testing::FaultKinSetup::_configureFaultFields(int spaceDim) {
    for (size_t iFault = 0; iFault < _faultObjects.size(); ++iFault) {
        // Set auxFieldDB on each kinematic source
        for (auto* src : _kinSrcs[iFault]) {
            src->auxFieldDB(&_faultAuxDB);
        }

        // Set discretization on fault auxiliary subfields
        for (size_t j = 0; j < _faultAuxSubfields.size(); ++j) {
            const auto& disc = _faultAuxDiscretizations[j];
            _faultObjects[iFault]->setAuxiliarySubfieldDiscretization(
                _faultAuxSubfields[j].c_str(),
                disc.basisOrder, disc.quadOrder, spaceDim - 1,
                disc.cellBasis, disc.feSpace, disc.isBasisContinuous);
        }
    }
}

void
pylith::testing::FaultKinSetup::_addLagrangeSubfields(
    pylith::problems::SolutionFactory& factory) {
    for (const auto& disc : _lagrangeDiscretizations) {
        factory.addLagrangeMultiplierFault(disc);
    }
}

std::vector<pylith::faults::FaultCohesive*>
pylith::testing::FaultKinSetup::_getFaultKinPtrs(void) {
    std::vector<pylith::faults::FaultCohesive*> ptrs;
    ptrs.reserve(_faultObjects.size());
    for (auto* f : _faultObjects) { ptrs.push_back(f); }
    return ptrs;
}
```

#### `_registerFaultExactSolution`

This replicates the logic in the current `TestFaultKin::_setExactSolution()`:

```cpp
void
pylith::testing::FaultKinSetup::_registerFaultExactSolution(
    PetscDM dm,
    const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnFns,
    const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnDotFns,
    size_t numDomainSubfields) {

    REQUIRE(!_faultObjects.empty());

    // Get cohesive cell DS
    PetscDMLabel label = nullptr;
    PetscIS is = nullptr;
    PetscInt cohesiveCell = -1;
    PylithCallPetsc(DMGetLabel(dm, pylith::topology::Mesh::cells_label_name, &label));
    PylithCallPetsc(DMLabelGetStratumIS(
        label, _faultObjects[0]->getCohesiveLabelValue(), &is));
    PylithCallPetsc(ISGetMinMax(is, &cohesiveCell, nullptr));
    PylithCallPetsc(ISDestroy(&is));

    PetscDS ds = nullptr;
    PylithCallPetsc(DMGetCellDS(dm, cohesiveCell, &ds, nullptr));

    // Register domain subfields on cohesive DS (without DM context)
    for (size_t i = 0; i < numDomainSubfields; ++i) {
        PylithCallPetsc(PetscDSSetExactSolution(ds, i, allExactSolnFns[i], nullptr));
        if (i < allExactSolnDotFns.size() && allExactSolnDotFns[i]) {
            PylithCallPetsc(PetscDSSetExactSolutionTimeDerivative(
                ds, i, allExactSolnDotFns[i], nullptr));
        }
    }

    // Register fault subfields (Lagrange multiplier)
    for (size_t i = 0; i < _numFaultSolnSubfields; ++i) {
        const size_t idx = numDomainSubfields + i;
        PylithCallPetsc(PetscDSSetExactSolution(
            ds, idx, allExactSolnFns[idx], nullptr));
    }
}
```

### 6.4 Multiple Kinematic Sources (Future Extension)

The `_kinSrcs` member is `std::vector<std::vector<KinSrc*>>`. Currently:
- Outer vector: one entry per fault.
- Inner vector: one entry (single `KinSrcStep`).

To support multiple sources per fault in the future:
- Add more entries to the inner vector.
- The rupture names array in `_transformFaultTopology` grows accordingly.
- `setEqRuptures` already accepts arrays, so no API change is needed.
- Each source can have a different `auxFieldDB` if needed (currently they share
  `_faultAuxDB`; this can be extended to per-source DBs).

---

## 7. Concrete Test Case Pattern

### 7.1 Structure

Each concrete test case class:
1. Inherits from `MMSTestCase` + appropriate mixin(s).
2. Defines analytical functions as **private static methods**.
3. Configures all parameters in its **constructor**.
4. Overrides the pure virtual hooks by delegating to mixin helpers.
5. Exposes **static factory methods** per mesh/discretization variant.

### 7.2 Example: No-Fault Linear Elasticity (`UniformStrain2D`)

**Header:** `tests/mmstests/linearelasticity/nofaults-2d/UniformStrain2D.hh`

```cpp
#pragma once

#include "tests/src/MMSTestCase.hh"
#include "tests/src/ElasticitySetup.hh"

namespace pylith {
    class UniformStrain2D;
}

class pylith::UniformStrain2D :
    public pylith::testing::MMSTestCase,
    public pylith::testing::ElasticitySetup {
public:

    UniformStrain2D(void);
    ~UniformStrain2D(void) = default;

    // Variant factory methods (return by value)
    static UniformStrain2D TriP1(void);
    static UniformStrain2D TriP2(void);
    static UniformStrain2D TriP3(void);
    static UniformStrain2D QuadQ1(void);
    static UniformStrain2D QuadQ2(void);
    static UniformStrain2D QuadQ3(void);

protected:

    // Hook overrides — delegate to mixin
    void _setupMaterials(void) override;
    void _setupSolutionSubfields(pylith::problems::SolutionFactory& factory) override;
    std::vector<pylith::materials::Material*> _getMaterialPtrs(void) override;

private:

    // Analytical functions
    static double density(double x, double y);
    static double vs(double x, double y);
    static double vp(double x, double y);
    static double disp_x(double x, double y);
    static double disp_y(double x, double y);

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim, PetscReal t,
                                          const PetscReal x[], PetscInt Nc,
                                          PetscScalar* s, void* ctx);
}; // class UniformStrain2D
```

**Implementation:** `tests/mmstests/linearelasticity/nofaults-2d/UniformStrain2D.cc`

```cpp
#include "UniformStrain2D.hh"
#include "pylith/problems/SolutionFactory.hh"
#include "pylith/scales/ElasticityScales.hh"
#include "pylith/bc/DirichletUserFn.hh"
#include "catch2/catch_test_macros.hpp"

// --- Constants ---
static const double STRAIN_XX = 0.1;
static const double STRAIN_YY = 0.25;
static const double STRAIN_XY = 0.3;

// --- Analytical functions ---
double pylith::UniformStrain2D::density(double x, double y) { return 2500.0; }
double pylith::UniformStrain2D::vs(double x, double y) { return 3000.0; }
double pylith::UniformStrain2D::vp(double x, double y) {
    return sqrt(3.0) * vs(x, y);
}

double pylith::UniformStrain2D::disp_x(double x, double y) {
    return STRAIN_XX * x + STRAIN_XY * y;
}
double pylith::UniformStrain2D::disp_y(double x, double y) {
    return STRAIN_XY * x + STRAIN_YY * y;
}

PetscErrorCode
pylith::UniformStrain2D::solnkernel_disp(PetscInt spaceDim, PetscReal t,
                                          const PetscReal x[], PetscInt Nc,
                                          PetscScalar* s, void* ctx) {
    s[0] = disp_x(x[0], x[1]);
    s[1] = disp_y(x[0], x[1]);
    return PETSC_SUCCESS;
}

// --- Constructor ---
pylith::UniformStrain2D::UniformStrain2D(void) {
    GenericComponent::setName("UniformStrain2D");

    // MMSTestCase members
    _spaceDim = 2;
    _cs.setSpaceDim(_spaceDim);
    _formulation = pylith::problems::Physics::QUASISTATIC;
    _isJacobianLinear = true;
    _enabledChecks = CHECK_ALL;

    const double lengthScale = 8.0e+3;
    pylith::scales::ElasticityScales::setQuasistaticElasticity(&_scales, lengthScale);

    // ElasticitySetup members
    _matAuxDB.addValue("density", density, "kg/m**3");
    _matAuxDB.addValue("vp", vp, "m/s");
    _matAuxDB.addValue("vs", vs, "m/s");
    _matAuxDB.setCoordSys(_cs);

    _matAuxSubfields = {"density", "shear_modulus", "bulk_modulus"};
    _matAuxDiscretizations = {
        pylith::topology::Field::Discretization(0, 1),  // density
        pylith::topology::Field::Discretization(0, 1),  // shear_modulus
        pylith::topology::Field::Discretization(0, 1),  // bulk_modulus
    };

    _matConfigs = {{"material-id=24", 24, false, false}};
    _rheology.useReferenceState(false);

    // Boundary conditions
    static const PylithInt constrainedDOF[2] = {0, 1};
    auto* bc = new pylith::bc::DirichletUserFn();
    bc->setSubfieldName("displacement");
    bc->setLabelName("boundary");
    bc->setLabelValue(1);
    bc->setConstrainedDOF(constrainedDOF, 2);
    bc->setUserFn(solnkernel_disp);
    _bcs.push_back(bc);

    // Exact solution
    _exactSolnFns = {solnkernel_disp};
}

// --- Hook overrides ---
void pylith::UniformStrain2D::_setupMaterials(void) {
    _createElasticityMaterials(_spaceDim, _formulation);
}

void pylith::UniformStrain2D::_setupSolutionSubfields(
    pylith::problems::SolutionFactory& factory) {
    factory.addDisplacement(_solnDiscretizations[0]);
    // No velocity for quasistatic
}

std::vector<pylith::materials::Material*>
pylith::UniformStrain2D::_getMaterialPtrs(void) {
    return _getElasticityMaterialPtrs();
}

// --- Variant factories ---
pylith::UniformStrain2D
pylith::UniformStrain2D::TriP1(void) {
    UniformStrain2D test;
    test._meshFilename = "data/tri.msh";
    test._useAsciiMesh = false;
    test._solnDiscretizations = {
        pylith::topology::Field::Discretization(1, 1),  // displacement
    };
    return test;
}

pylith::UniformStrain2D
pylith::UniformStrain2D::TriP2(void) {
    UniformStrain2D test;
    test._meshFilename = "data/tri.msh";
    test._useAsciiMesh = false;
    test._solnDiscretizations = {
        pylith::topology::Field::Discretization(2, 2),
    };
    return test;
}

// QuadQ1, QuadQ2, QuadQ3, TriP3 follow the same pattern...
```

### 7.3 Example: Fault Test Case (`TwoBlocksStatic`)

**Header:** `tests/mmstests/linearelasticity/faults-2d/TwoBlocksStatic.hh`

```cpp
#pragma once

#include "tests/src/MMSTestCase.hh"
#include "tests/src/ElasticitySetup.hh"
#include "tests/src/FaultKinSetup.hh"

namespace pylith {
    class TwoBlocksStatic;
}

class pylith::TwoBlocksStatic :
    public pylith::testing::MMSTestCase,
    public pylith::testing::ElasticitySetup,
    public pylith::testing::FaultKinSetup {
public:
    TwoBlocksStatic(void);
    ~TwoBlocksStatic(void) = default;

    static TwoBlocksStatic TriP1(void);
    static TwoBlocksStatic TriP2(void);
    static TwoBlocksStatic TriP3(void);
    static TwoBlocksStatic QuadQ1(void);
    static TwoBlocksStatic QuadQ2(void);
    static TwoBlocksStatic QuadQ3(void);

protected:
    void _setupFaultTopology(void) override;
    void _setupMaterials(void) override;
    void _setupFaultFields(void) override;
    void _setupSolutionSubfields(pylith::problems::SolutionFactory& factory) override;
    std::vector<pylith::materials::Material*> _getMaterialPtrs(void) override;
    std::vector<pylith::faults::FaultCohesive*> _getFaultPtrs(void) override;
    void _registerExactSolution(PetscDS ds, PetscDM dm) override;
    void _registerExactSolutionFault(PetscDM dm) override;

private:
    static const double LENGTH_SCALE;
    static const double TIME_SCALE;
    static const double RIGIDITY_SCALE;
    static const double AMPLITUDE;
    static const double X_FAULT;

    static double density(double x, double y);
    static double vs(double x, double y);
    static double vp(double x, double y);
    static double initiation_time(double x, double y);
    static double finalslip_opening(double x, double y);
    static double finalslip_leftlateral(double x, double y);
    static double disp_x(double x, double y);
    static double disp_y(double x, double y, PetscInt flag);
    static double faulttraction_x(double x, double y);
    static double faulttraction_y(double x, double y);

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim, PetscReal t,
                                          const PetscReal x[], PetscInt Nc,
                                          PetscScalar* s, void* ctx);
    static PetscErrorCode solnkernel_lagrangemultiplier(PetscInt spaceDim,
                                                        PetscReal t,
                                                        const PetscReal x[],
                                                        PetscInt Nc,
                                                        PetscScalar* s, void* ctx);
};
```

**Implementation key parts:**

```cpp
pylith::TwoBlocksStatic::TwoBlocksStatic(void) {
    GenericComponent::setName("TwoBlocksStatic");

    _spaceDim = 2;
    _cs.setSpaceDim(_spaceDim);
    _formulation = pylith::problems::Physics::QUASISTATIC;
    _isJacobianLinear = true;

    _scales.setLengthScale(LENGTH_SCALE);
    _scales.setTimeScale(TIME_SCALE);
    _scales.setRigidityScale(RIGIDITY_SCALE);

    // Material config (3 blocks)
    _matConfigs = {
        {"material-id=10", 10, false, false},
        {"material-id=20", 20, false, false},
        {"material-id=15", 15, false, false},
    };
    _matAuxDB.addValue("density", density, "kg/m**3");
    _matAuxDB.addValue("vp", vp, "m/s");
    _matAuxDB.addValue("vs", vs, "m/s");
    _matAuxDB.setCoordSys(_cs);

    _matAuxSubfields = {"density", "shear_modulus", "bulk_modulus"};
    _matAuxDiscretizations = {
        {0, 1}, {0, 1}, {0, 1},
    };

    // Fault config
    _faultConfigs = {{100, "fault_xpos_faces", 1, "", 0}};
    _numFaultSolnSubfields = 1;

    // Kinematic source
    auto* kinSrc = new pylith::faults::KinSrcStep();
    kinSrc->setOriginTime(0.0);
    _kinSrcs = {{kinSrc}};

    _faultAuxDB.addValue("initiation_time", initiation_time, "s");
    _faultAuxDB.addValue("final_slip_opening", finalslip_opening, "m");
    _faultAuxDB.addValue("final_slip_left_lateral", finalslip_leftlateral, "m");
    _faultAuxDB.setCoordSys(_cs);

    _faultAuxSubfields = {"slip"};
    _faultAuxDiscretizations = {{0, 1}};

    // Boundary conditions
    static const PylithInt dof[2] = {0, 1};
    {
        auto* bc = new pylith::bc::DirichletUserFn();
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary_xpos");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(dof, 2);
        bc->setUserFn(solnkernel_disp);
        _bcs.push_back(bc);
    }
    {
        auto* bc = new pylith::bc::DirichletUserFn();
        bc->setSubfieldName("displacement");
        bc->setLabelName("boundary_xneg");
        bc->setLabelValue(1);
        bc->setConstrainedDOF(dof, 2);
        bc->setUserFn(solnkernel_disp);
        _bcs.push_back(bc);
    }

    // Exact solution (domain + fault)
    _exactSolnFns = {solnkernel_disp, solnkernel_lagrangemultiplier};
}

// Hook overrides
void pylith::TwoBlocksStatic::_setupFaultTopology(void) {
    _transformFaultTopology(_mesh);
}
void pylith::TwoBlocksStatic::_setupMaterials(void) {
    _createElasticityMaterials(_spaceDim, _formulation);
}
void pylith::TwoBlocksStatic::_setupFaultFields(void) {
    _configureFaultFields(_spaceDim);
}
void pylith::TwoBlocksStatic::_setupSolutionSubfields(
    pylith::problems::SolutionFactory& factory) {
    factory.addDisplacement(_solnDiscretizations[0]);
    _addLagrangeSubfields(factory);
}
std::vector<pylith::materials::Material*>
pylith::TwoBlocksStatic::_getMaterialPtrs(void) {
    return _getElasticityMaterialPtrs();
}
std::vector<pylith::faults::FaultCohesive*>
pylith::TwoBlocksStatic::_getFaultPtrs(void) {
    return _getFaultKinPtrs();
}
void pylith::TwoBlocksStatic::_registerExactSolution(PetscDS ds, PetscDM dm) {
    // Domain subfields only (with DM context for displacement kernel)
    const size_t numDomain = _exactSolnFns.size() - _numFaultSolnSubfields;
    for (size_t i = 0; i < numDomain; ++i) {
        PylithCallPetsc(PetscDSSetExactSolution(ds, i, _exactSolnFns[i], dm));
    }
}
void pylith::TwoBlocksStatic::_registerExactSolutionFault(PetscDM dm) {
    const size_t numDomain = _exactSolnFns.size() - _numFaultSolnSubfields;
    _registerFaultExactSolution(dm, _exactSolnFns, _exactSolnDotFns, numDomain);
}
```

### 7.4 Example: Poroelasticity (`PressureGradient`)

```cpp
class pylith::PressureGradient :
    public pylith::testing::MMSTestCase,
    public pylith::testing::PoroelasticitySetup {
public:
    PressureGradient(void);

    static PressureGradient TriP2(void);
    // ...

protected:
    void _setupMaterials(void) override {
        _createPoroelasticityMaterials(_spaceDim, _formulation);
    }
    void _setupSolutionSubfields(pylith::problems::SolutionFactory& factory) override {
        factory.addDisplacement(_solnDiscretizations[0]);
        factory.addPressure(_solnDiscretizations[1]);
        factory.addTraceStrain(_solnDiscretizations[2]);
        // For time-dependent with 6 subfields:
        // factory.addVelocity(_solnDiscretizations[3]);
        // factory.addPressureDot(_solnDiscretizations[4]);
        // factory.addTraceStrainDot(_solnDiscretizations[5]);
    }
    std::vector<pylith::materials::Material*> _getMaterialPtrs(void) override {
        return _getPoroelasticityMaterialPtrs();
    }

private:
    // Analytical functions for displacement, pressure, trace strain,
    // material properties (porosity, permeability, fluid viscosity, etc.)
    // ...
};
```

### 7.5 Dynamic Formulation Example (`PlanePWave2D`)

```cpp
pylith::PlanePWave2D::PlanePWave2D(void) {
    // ...
    _formulation = pylith::problems::Physics::DYNAMIC;
    _enabledChecks = CHECK_DISCRETIZATION | CHECK_RESIDUAL;  // No Jacobian checks
    // ...
    _exactSolnFns = {solnkernel_disp, solnkernel_vel};
    _exactSolnDotFns = {solnkernel_disp_dot, solnkernel_vel_dot};
}

void pylith::PlanePWave2D::_setupSolutionSubfields(
    pylith::problems::SolutionFactory& factory) {
    factory.addDisplacement(_solnDiscretizations[0]);
    factory.addVelocity(_solnDiscretizations[1]);
}
```

---

## 8. `TestCases.cc` and the `MMS_TEST_CASE` Macro

### 8.1 Macro Definition

**File:** `tests/src/MMSTestMacros.hh`

```cpp
#pragma once

#include "catch2/catch_test_macros.hpp"

/// Register all four MMS checks for a test case variant.
/// Checks are skipped at runtime if disabled via _enabledChecks.
#define MMS_TEST_CASE(TestClass, Variant)                                     \
    TEST_CASE(#TestClass "::" #Variant "::testDiscretization",                \
              "[" #TestClass "][" #Variant "][discretization]") {              \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_DISCRETIZATION))           \
            test.testDiscretization();                                         \
    }                                                                         \
    TEST_CASE(#TestClass "::" #Variant "::testResidual",                      \
              "[" #TestClass "][" #Variant "][residual]") {                    \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_RESIDUAL))                 \
            test.testResidual();                                               \
    }                                                                         \
    TEST_CASE(#TestClass "::" #Variant "::testJacobianTaylorSeries",          \
              "[" #TestClass "][" #Variant "][Jacobian Taylor series]") {      \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_JACOBIAN_TAYLOR))          \
            test.testJacobianTaylorSeries();                                   \
    }                                                                         \
    TEST_CASE(#TestClass "::" #Variant "::testJacobianFiniteDiff",            \
              "[" #TestClass "][" #Variant "][Jacobian finite difference]") {  \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_JACOBIAN_FD))              \
            test.testJacobianFiniteDiff();                                     \
    }
```

### 8.2 Usage in `TestCases.cc`

```cpp
// tests/mmstests/linearelasticity/nofaults-2d/TestCases.cc

#include "tests/src/MMSTestMacros.hh"
#include "UniformStrain2D.hh"
#include "Gravity2D.hh"
#include "GravityRefState2D.hh"
#include "BodyForce2D.hh"
#include "RigidBodyAcc2D.hh"
#include "PlanePWave2D.hh"

MMS_TEST_CASE(UniformStrain2D, TriP1)
MMS_TEST_CASE(UniformStrain2D, TriP2)
MMS_TEST_CASE(UniformStrain2D, TriP3)
MMS_TEST_CASE(UniformStrain2D, QuadQ1)
MMS_TEST_CASE(UniformStrain2D, QuadQ2)
MMS_TEST_CASE(UniformStrain2D, QuadQ3)

MMS_TEST_CASE(Gravity2D, TriP2)
MMS_TEST_CASE(Gravity2D, TriP3)
MMS_TEST_CASE(Gravity2D, QuadQ2)
MMS_TEST_CASE(Gravity2D, QuadQ3)

MMS_TEST_CASE(GravityRefState2D, TriP1)
MMS_TEST_CASE(GravityRefState2D, TriP2)
MMS_TEST_CASE(GravityRefState2D, TriP3)
MMS_TEST_CASE(GravityRefState2D, QuadQ1)
MMS_TEST_CASE(GravityRefState2D, QuadQ2)
MMS_TEST_CASE(GravityRefState2D, QuadQ3)

// ... etc.
```

The `GravityRefState2D` constructor sets
`_enabledChecks = CHECK_DISCRETIZATION | CHECK_RESIDUAL | CHECK_JACOBIAN_TAYLOR`
(omitting `CHECK_JACOBIAN_FD`), so the finite-difference test is skipped at runtime.

---

## 9. File Layout

### 9.1 New Shared Infrastructure (`tests/src/`)

```
tests/src/
    MMSTest.hh                       (UNCHANGED)
    MMSTest.cc                       (UNCHANGED)
    MMSTestCase.hh                   (NEW)
    MMSTestCase.cc                   (NEW)
    MMSTestMacros.hh                 (NEW)
    ElasticitySetup.hh               (NEW)
    ElasticitySetup.cc               (NEW)
    IncompressibleElasticitySetup.hh (NEW)
    IncompressibleElasticitySetup.cc (NEW)
    PoroelasticitySetup.hh           (NEW)
    PoroelasticitySetup.cc           (NEW)
    FaultKinSetup.hh                 (NEW)
    FaultKinSetup.cc                 (NEW)
    driver_catch2.cc                 (UNCHANGED)
    ... (other existing files unchanged)
```

### 9.2 Per-Directory Structure (e.g., `linearelasticity/nofaults-2d/`)

```
tests/mmstests/linearelasticity/nofaults-2d/
    UniformStrain2D.hh      (REWRITTEN)
    UniformStrain2D.cc      (REWRITTEN)
    Gravity2D.hh            (REWRITTEN)
    Gravity2D.cc            (REWRITTEN)
    GravityRefState2D.hh    (REWRITTEN)
    GravityRefState2D.cc    (REWRITTEN)
    BodyForce2D.hh          (REWRITTEN)
    BodyForce2D.cc          (REWRITTEN)
    RigidBodyAcc2D.hh       (REWRITTEN)
    RigidBodyAcc2D.cc       (REWRITTEN)
    PlanePWave2D.hh         (REWRITTEN)
    PlanePWave2D.cc         (REWRITTEN)
    TestCases.cc            (SIMPLIFIED using macro)
    Makefile.am             (UPDATED — removes TestLinearElasticity sources)
    data/                   (UNCHANGED)
```

**Removed files:**
- `TestLinearElasticity.hh` / `.cc` — replaced by `MMSTestCase` + `ElasticitySetup`
- `TestFaultKin.hh` / `.cc` — replaced by `MMSTestCase` + `ElasticitySetup` + `FaultKinSetup`
- `TestIncompressibleElasticity.hh` / `.cc` — replaced by `MMSTestCase` + `IncompressibleElasticitySetup`
- `TestLinearPoroelasticity.hh` / `.cc` — replaced by `MMSTestCase` + `PoroelasticitySetup`

### 9.3 Updated `Makefile.am` Pattern

```makefile
include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_linearelasticity_nofaults2d

TESTS = \
    run_UniformStrain2D.sh \
    run_Gravity2D.sh \
    run_GravityRefState2D.sh \
    run_BodyForce2D.sh \
    run_RigidBodyAcc2D.sh \
    run_PlanePWave2D.sh

check_SCRIPTS = $(TESTS)
check_PROGRAMS = $(MMS_DRIVER)

mmstest_linearelasticity_nofaults2d_SOURCES = \
    $(top_srcdir)/tests/src/MMSTest.cc \
    $(top_srcdir)/tests/src/MMSTestCase.cc \
    $(top_srcdir)/tests/src/ElasticitySetup.cc \
    $(top_srcdir)/tests/src/driver_catch2.cc \
    TestCases.cc \
    UniformStrain2D.cc \
    Gravity2D.cc \
    GravityRefState2D.cc \
    BodyForce2D.cc \
    RigidBodyAcc2D.cc \
    PlanePWave2D.cc

# For fault directories, add FaultKinSetup.cc:
# $(top_srcdir)/tests/src/FaultKinSetup.cc \

run_%.sh:
    echo "#!/bin/bash" > $@
    echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
    chmod +x $@

dist_noinst_HEADERS = \
    UniformStrain2D.hh \
    Gravity2D.hh \
    GravityRefState2D.hh \
    BodyForce2D.hh \
    RigidBodyAcc2D.hh \
    PlanePWave2D.hh
```

---

## 10. Initialization Sequence Diagram

```
TestCases.cc:  UniformStrain2D::TriP2().testResidual()
    │
    ▼
MMSTest::testResidual()        [unchanged from current code]
    │
    ├── _initialize()          [dispatches to MMSTestCase::_initialize()]
    │       │
    │       ├── _readMesh()                  [MMSTestCase default: MeshIOPetsc]
    │       │
    │       ├── _setupFaultTopology()        [no-op, no FaultKinSetup]
    │       │
    │       ├── _nondimensionalize()         [MMSTestCase default]
    │       │
    │       ├── _setupMaterials()            [UniformStrain2D delegates to
    │       │                                 ElasticitySetup::_createElasticityMaterials()]
    │       │
    │       ├── _setupFaultFields()          [no-op, no FaultKinSetup]
    │       │
    │       ├── _setupProblem()              [MMSTestCase default:
    │       │       │                         calls _getMaterialPtrs(), _getFaultPtrs()]
    │       │       ├── _problem->setMaterials(...)
    │       │       ├── _problem->setBoundaryConditions(...)
    │       │       ├── _problem->setStartTime/EndTime/InitialTimeStep/Formulation
    │       │       └── (no faults to set)
    │       │
    │       ├── Create _solution, SolutionFactory
    │       ├── _setupSolutionSubfields()    [UniformStrain2D:
    │       │                                 factory.addDisplacement(...)]
    │       ├── _problem->setSolution(...)
    │       │
    │       └── MMSTest::_initialize()       [PETSc TS/SNES setup, calls _setExactSolution()]
    │               │
    │               └── _setExactSolution()  [MMSTestCase:
    │                       │                 DMGetDS, _registerExactSolution(),
    │                       │                 _registerExactSolutionFault()]
    │                       ├── PetscDSSetExactSolution(ds, 0, solnkernel_disp, nullptr)
    │                       └── (no fault registration)
    │
    └── DMTSCheckResidual(...)               [PETSc MMS verification]
```

For a fault test (`TwoBlocksStatic`), the diagram adds:

```
    ├── _setupFaultTopology()        [TwoBlocksStatic delegates to
    │                                 FaultKinSetup::_transformFaultTopology()]
    │       └── For each fault: create FaultCohesiveKin, setEqRuptures,
    │           transformTopology (inserts cohesive cells)
    ...
    ├── _setupFaultFields()          [TwoBlocksStatic delegates to
    │                                 FaultKinSetup::_configureFaultFields()]
    │       └── For each fault/kinSrc: set auxFieldDB, set discretizations
    ...
    ├── _setupProblem()
    │       ├── _problem->setMaterials(...)
    │       ├── _problem->setInterfaces(faults...)    ← from _getFaultPtrs()
    │       └── _problem->setBoundaryConditions(...)
    ...
    ├── _setupSolutionSubfields()
    │       ├── factory.addDisplacement(...)
    │       └── factory.addLagrangeMultiplierFault(...)  ← from _addLagrangeSubfields()
    ...
    └── _setExactSolution()
            ├── _registerExactSolution(ds, dm)
            │       └── PetscDSSetExactSolution(ds, 0, solnkernel_disp, dm)
            │           (note: passes DM as context for cell-side detection)
            └── _registerExactSolutionFault(dm)
                    └── FaultKinSetup::_registerFaultExactSolution(...)
                            ├── Get cohesive cell DS
                            ├── Register domain fns without DM context
                            └── Register Lagrange multiplier fn
```

---

## 11. Design Decisions and Rationale

### 11.1 Why Multiple Inheritance Instead of Composition

The equation setup and fault setup are **orthogonal concerns** that naturally
compose. Multiple inheritance maps directly to the "is-a" relationship:
- `TwoBlocksStatic` **is** an MMS test case, **uses** elasticity setup, **uses** fault setup.

Composition (having a pointer to an equation helper) would require forwarding
calls and passing `this` context, adding boilerplate without benefit.

### 11.2 Why Virtual Hooks Instead of Template Methods (CRTP)

- The project uses C-style C++ with minimal template usage.
- Virtual dispatch is familiar to contributors.
- The cost of virtual calls is negligible (called once during test initialization).
- No template instantiation compile-time overhead.

### 11.3 Why Test Cases Return by Value

The variant factories (e.g., `TriP1()`) return the test object **by value**.
This is safe because:
- `MMSTestCase` and mixins own heap-allocated objects (BCs, materials, faults)
  via `std::vector` of pointers. Move semantics handle transfer.
- The test object is constructed, used for a single test method call, then
  destroyed. No sharing or copying occurs.

If move semantics are problematic (due to raw pointer ownership), an alternative
is returning `std::unique_ptr<TestClass>`. However, value semantics are simpler
and match the existing pattern where `createData()` returns a raw pointer to
heap-allocated data.

### 11.4 Why `_enabledChecks` Bitmask Instead of Separate Booleans

A bitmask is compact, extensible (add new check types without adding members),
and allows the macro to query any check uniformly.

### 11.5 Diamond Inheritance Prevention

`MMSTestCase` does not inherit from any mixin. Equation mixins and
`FaultKinSetup` do not inherit from `MMSTestCase`. There is no diamond:

```
MMSTest ← MMSTestCase ← ConcreteTest
                              ↑ (also inherits)
              ElasticitySetup ┘
              FaultKinSetup   ┘
```

The mixins are **non-polymorphic helper classes** — they have no virtual base
and no relationship to `MMSTest` or `MMSTestCase` in the class hierarchy.
The concrete test class is the only class with multiple base classes, and those
bases are unrelated to each other.

### 11.6 Ownership Model

| Object | Owner | Lifetime |
|--------|-------|----------|
| `Mesh` | `MMSTest` (base) | Entire test |
| `TimeDependent` problem | `MMSTest` (base) | Entire test |
| `Solution` field | `MMSTest` (base) | Entire test |
| `BoundaryCondition*` objects | `MMSTestCase::_bcs` | Destroyed in `~MMSTestCase` |
| `GravityField*` | `MMSTestCase::_gravityField` | Destroyed in `~MMSTestCase` |
| `Elasticity*` materials | `ElasticitySetup::_elasticityMaterials` | Destroyed in `~ElasticitySetup` |
| `IsotropicLinearElasticity` rheology | `ElasticitySetup::_rheology` | Stack, lives with mixin |
| `FaultCohesiveKin*` objects | `FaultKinSetup::_faultObjects` | Destroyed in `~FaultKinSetup` |
| `KinSrc*` objects | `FaultKinSetup::_kinSrcs` | Destroyed in `~FaultKinSetup` |
| `UserFunctionDB` objects | Stack members in mixin | Lives with mixin |

Destruction order: C++ destroys members in reverse declaration order, and base
classes in reverse inheritance order. Since the problem holds raw pointers to
materials/faults/BCs (not owning), the test infrastructure must ensure mixins
outlive the `MMSTest` base. This is guaranteed because `MMSTestCase` (which
inherits from `MMSTest`) is destroyed before the mixin bases in reverse order:
`~FaultKinSetup` → `~ElasticitySetup` → `~MMSTestCase` → `~MMSTest`.

Wait — C++ destroys in the order: most-derived destructor first, then base
classes in reverse declaration order. For
`class Concrete : public MMSTestCase, public ElasticitySetup, public FaultKinSetup`,
destruction order is:
1. `~Concrete()`
2. `~FaultKinSetup()`
3. `~ElasticitySetup()`
4. `~MMSTestCase()`
5. `~MMSTest()`

This means `MMSTest::~MMSTest()` (which calls `VecDestroy` on exact solution
vectors) runs **after** the mixin destructors free materials/faults. This is
safe because `MMSTest::~MMSTest()` only destroys PETSc vectors it owns directly,
not the materials/faults.

The `TimeDependent` problem is destroyed in `MMSTest::~MMSTest()`. The problem
holds raw non-owning pointers to materials and faults. By the time `~MMSTest`
runs, those objects are already freed. However, `~TimeDependent` does NOT delete
the materials/faults it points to (it stores non-owning pointers set via
`setMaterials`/`setInterfaces`). So there is no double-free or use-after-free.

---

## 12. Migration Plan

### 12.1 Phase 1: Create Shared Infrastructure

1. Implement `MMSTestCase.hh` / `.cc` in `tests/src/`.
2. Implement `ElasticitySetup.hh` / `.cc` in `tests/src/`.
3. Implement `FaultKinSetup.hh` / `.cc` in `tests/src/`.
4. Implement `MMSTestMacros.hh` in `tests/src/`.
5. Update `tests/src/Makefile.am` to include new files.

### 12.2 Phase 2: Migrate Linear Elasticity No-Fault 2D

1. Rewrite `UniformStrain2D.hh` / `.cc` using new pattern.
2. Update `TestCases.cc` to use `MMS_TEST_CASE` macro.
3. Remove `TestLinearElasticity.hh` / `.cc`.
4. Update `Makefile.am`.
5. Verify all tests pass.

### 12.3 Phase 3: Migrate Remaining Directories

Proceed one directory at a time:
1. `linearelasticity/faults-2d` — requires `FaultKinSetup`.
2. `linearelasticity/nofaults-3d` — straightforward, same as 2D but `_spaceDim=3`.
3. `incompressibleelasticity/nofaults-2d` — requires `IncompressibleElasticitySetup`.
4. `poroelasticity/nofaults-2d` — requires `PoroelasticitySetup`.

### 12.4 Phase 4: Add New Equation Mixins

Implement `IncompressibleElasticitySetup` and `PoroelasticitySetup` as each
directory is migrated. These are straightforward copies of the `ElasticitySetup`
pattern with different material/rheology types.

### 12.5 Validation

After each phase, run the full test suite to confirm:
- All existing MMS tests still pass.
- Test names and Catch2 tags are preserved (for CI filtering).
- No regressions in test coverage.

---

## 13. API Reference Summary

### 13.1 Key External APIs Used

| Class | Method | Purpose |
|-------|--------|---------|
| `SolutionFactory` | `addDisplacement(Discretization)` | Register displacement subfield |
| `SolutionFactory` | `addVelocity(Discretization)` | Register velocity subfield |
| `SolutionFactory` | `addPressure(Discretization)` | Register pressure subfield |
| `SolutionFactory` | `addPressureDot(Discretization)` | Register pressure time derivative |
| `SolutionFactory` | `addTraceStrain(Discretization)` | Register trace strain subfield |
| `SolutionFactory` | `addTraceStrainDot(Discretization)` | Register trace strain time derivative |
| `SolutionFactory` | `addTemperature(Discretization)` | Register temperature subfield |
| `SolutionFactory` | `addLagrangeMultiplierFault(Discretization)` | Register fault Lagrange multiplier |
| `DirichletUserFn` | `setUserFn(PetscUserFieldFunc)` | Set analytical BC function |
| `DirichletUserFn` | `setUserFnDot(PetscUserFieldFunc)` | Set analytical BC time derivative |
| `DirichletUserFn` | `setConstrainedDOF(int*, int)` | Set constrained DOFs |
| `DirichletUserFn` | `setSubfieldName(const char*)` | Target solution subfield |
| `DirichletUserFn` | `setLabelName(const char*)` | Mesh label for BC region |
| `DirichletUserFn` | `setLabelValue(int)` | Label value for BC region |
| `NeumannUserFn` | `setUserFn(PetscBdPointFn*)` | Set boundary integration kernel |
| `NeumannUserFn` | `setSubfieldName(const char*)` | Target solution subfield |
| `NeumannUserFn` | `setLabelName(const char*)` | Mesh label for BC region |
| `FaultCohesiveKin` | `setCohesiveLabelValue(int)` | Label value for cohesive cells |
| `FaultCohesiveKin` | `setSurfaceLabelName(const char*)` | Surface label from mesh |
| `FaultCohesiveKin` | `setEqRuptures(names, n, srcs, n)` | Set kinematic sources |
| `FaultCohesiveKin` | `transformTopology(Mesh*)` | Insert cohesive cells |
| `FaultCohesiveKin` | `setAuxiliarySubfieldDiscretization(...)` | Set fault aux discretization |
| `KinSrc` | `setOriginTime(PylithReal)` | Origin time for rupture |
| `KinSrc` | `auxFieldDB(SpatialDB*)` | Set spatial DB for slip params |
| `UserFunctionDB` | `addValue(name, fn, units)` | Register analytical function |
| `UserFunctionDB` | `setCoordSys(CSCart)` | Set coordinate system |
| `Material` | `setAuxiliaryFieldDB(SpatialDB*)` | Set material property DB |
| `Material` | `setAuxiliarySubfieldDiscretization(...)` | Set aux discretization |
| `Material` | `setFormulation(FormulationEnum)` | Quasistatic or dynamic |
| `Material` | `setBulkRheology(RheologyElasticity*)` | Set rheology |
| `Material` | `setLabelValue(int)` | Material label in mesh |
| `Material` | `setName(const char*)` | Human-readable name |
| `Material` | `setIdentifier(const char*)` | String identifier |

### 13.2 PETSc Function Signatures

```cpp
// Exact solution / Dirichlet BC kernel
typedef PetscErrorCode (*PetscUserFieldFunc)(
    PetscInt dim, PetscReal t, const PetscReal x[],
    PetscInt Nc, PetscScalar u[], void* ctx);

// Neumann BC boundary integration kernel
typedef void (*PetscBdPointFn)(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[],
    PetscScalar f[]);
```

### 13.3 `Discretization` Struct Fields

```cpp
struct Discretization {
    int basisOrder;       // Polynomial order of basis functions
    int quadOrder;        // Quadrature order
    int dimension;        // Topological dimension (-1 = use default)
    int cellBasis;        // Cell basis type (DEFAULT_BASIS = default)
    int feSpace;          // Finite element space (POLYNOMIAL_SPACE = default)
    bool isBasisContinuous; // Whether basis is continuous (true = default)
};
```

Constructor shorthand: `Discretization(basisOrder, quadOrder)` uses defaults
for remaining fields.

---

## 14. Checklist for Adding a New Equation Type

To add support for a new governing equation (e.g., thermoplasticity):

1. **Create mixin files:**
   - `tests/src/ThermoplasticitySetup.hh`
   - `tests/src/ThermoplasticitySetup.cc`

2. **Define in the mixin:**
   - `MaterialConfig` struct with equation-specific flags.
   - `_create*Materials()` method that allocates the correct material/rheology.
   - `_get*MaterialPtrs()` method returning `vector<Material*>`.
   - `UserFunctionDB` and aux subfield vectors for material properties.
   - Owned material and rheology objects.

3. **Create test case directory:**
   - `tests/mmstests/thermoplasticity/nofaults-2d/`
   - Add `data/` subdirectory with mesh files.

4. **Write concrete test case classes:**
   - Inherit from `MMSTestCase` + `ThermoplasticitySetup` (+ `FaultKinSetup` if needed).
   - Define analytical functions for material properties and solution.
   - Implement `_setupMaterials()`, `_setupSolutionSubfields()`, `_getMaterialPtrs()`.
   - Expose variant factory methods.

5. **Write `TestCases.cc`:**
   - Use `MMS_TEST_CASE` macro for each variant.

6. **Write `Makefile.am`:**
   - Include `$(top_srcdir)/tests/src/MMSTestCase.cc` and
     `$(top_srcdir)/tests/src/ThermoplasticitySetup.cc`.

7. **No changes needed to:**
   - `MMSTest.hh` / `.cc`
   - `MMSTestCase.hh` / `.cc`
   - Any other mixin
   - Any other test directory

---

## 15. Checklist for Adding a New Test Case to an Existing Equation

To add a new analytical solution test (e.g., `NonuniformStrain2D` for linear
elasticity):

1. **Create files:**
   - `tests/mmstests/linearelasticity/nofaults-2d/NonuniformStrain2D.hh`
   - `tests/mmstests/linearelasticity/nofaults-2d/NonuniformStrain2D.cc`

2. **Define the class:**
   - Inherit from `MMSTestCase` + `ElasticitySetup`.
   - Write analytical functions (density, vs, vp, displacement components).
   - Write the PETSc solution kernel function(s).
   - In constructor: set up `_matAuxDB`, `_matConfigs`, `_bcs`, `_exactSolnFns`.
   - Implement hook overrides (typically 3 one-liners delegating to mixin).
   - Write variant factories (`TriP1`, `TriP2`, etc.).

3. **Update `TestCases.cc`:**
   - Add `#include "NonuniformStrain2D.hh"`.
   - Add `MMS_TEST_CASE(NonuniformStrain2D, TriP2)` etc.

4. **Update `Makefile.am`:**
   - Add `NonuniformStrain2D.cc` to `_SOURCES`.
   - Add `NonuniformStrain2D.hh` to `dist_noinst_HEADERS`.
   - Add `run_NonuniformStrain2D.sh` to `TESTS`.

---

## 16. Potential Issues and Mitigations

### 16.1 Copy/Move Semantics of Test Objects

The variant factories return by value. The test classes own heap-allocated
objects (BCs, materials, faults) via `std::vector<T*>`. Default copy would
double-free these pointers.

**Mitigation:** Either:
- (A) Delete copy constructor/assignment, implement move constructor/assignment.
- (B) Use `std::unique_ptr<T>` in the vectors instead of raw pointers.
- (C) Return by `std::unique_ptr<TestClass>` from factories.

**Recommendation:** Option (A) — implement move semantics. The vectors transfer
ownership cleanly via `std::move`. This matches the existing pattern where data
structs were heap-allocated and returned by raw pointer.

### 16.2 Static Data in Analytical Functions

Some test cases use `static` variables for constants (e.g., `STRAIN_XX`).
With the new design these become file-scope constants or `static constexpr`
members. No functional change, but ensure they are not shared across
translation units in ways that could conflict.

### 16.3 Thread Safety

MMS tests run sequentially within a single Catch2 process. No thread safety
concerns.

### 16.4 Compile Time

Each test directory compiles into a single executable. The shared infrastructure
files (`MMSTestCase.cc`, `ElasticitySetup.cc`) are compiled once per executable.
This is equivalent to the current pattern where `MMSTest.cc` and
`TestLinearElasticity.cc` are compiled per executable. No compile time regression.

---

## 17. Complete Implementation Specification for `MMSTestCase`

Below is the complete, definitive interface. The implementation agent should
produce files matching this interface exactly.

### Required includes for `MMSTestCase.cc`:

```cpp
#include "tests/src/MMSTestCase.hh"
#include "pylith/problems/TimeDependent.hh"
#include "pylith/problems/SolutionFactory.hh"
#include "pylith/topology/Mesh.hh"
#include "pylith/topology/MeshOps.hh"
#include "pylith/topology/Field.hh"
#include "pylith/meshio/MeshIOAscii.hh"
#include "pylith/meshio/MeshIOPetsc.hh"
#include "pylith/exceptions/error.hh"
#include "pylith/utils/journals.hh"
#include "spatialdata/spatialdb/GravityField.hh"
#include "catch2/catch_test_macros.hpp"
#include "petscdm.h"
#include "petscds.h"
```

### Methods that MUST remain virtual (overridable by concrete classes):

| Method | Default behavior | When overridden |
|--------|-----------------|-----------------|
| `_readMesh()` | Reads mesh via MeshIOAscii or MeshIOPetsc | Rarely |
| `_setupFaultTopology()` | No-op | By `FaultKinSetup` via concrete class |
| `_nondimensionalize()` | Sets coord sys, calls MeshOps::nondimensionalize | Rarely |
| `_setupMaterials()` | **PURE VIRTUAL** | Always (by equation mixin) |
| `_setupSolutionSubfields()` | **PURE VIRTUAL** | Always (by equation mixin) |
| `_setupFaultFields()` | No-op | By `FaultKinSetup` via concrete class |
| `_setupProblem()` | Wires problem object | Rarely |
| `_registerExactSolution()` | Iterates `_exactSolnFns` | By fault tests (need DM context) |
| `_registerExactSolutionFault()` | No-op | By `FaultKinSetup` via concrete class |
| `_getMaterialPtrs()` | Returns empty vector | Always (by equation mixin) |
| `_getFaultPtrs()` | Returns empty vector | By `FaultKinSetup` via concrete class |

---

## 18. Existing Test Cases to Migrate

### 18.1 `linearelasticity/nofaults-2d` (6 test cases)

| Test Case | Formulation | Gravity | Body Force | Ref State | Notes |
|-----------|-------------|---------|------------|-----------|-------|
| `UniformStrain2D` | Quasistatic | No | No | No | Uniform strain field; all 4 checks |
| `Gravity2D` | Quasistatic | Yes | No | No | Quadratic displacement; P2+ only |
| `GravityRefState2D` | Quasistatic | Yes | No | Yes | Zero displacement (ref state absorbs gravity); skip FD Jacobian |
| `BodyForce2D` | Quasistatic | No | Yes | No | Quadratic displacement; P2+ only |
| `RigidBodyAcc2D` | Dynamic | No | No | No | Rigid body acceleration; only discretization + residual checks |
| `PlanePWave2D` | Dynamic | No | No | No | Plane P-wave; only discretization + residual checks |

### 18.2 `linearelasticity/nofaults-3d` (4 test cases)

| Test Case | Formulation | Gravity | Body Force | Ref State |
|-----------|-------------|---------|------------|-----------|
| `UniformStrain3D` | Quasistatic | No | No | No |
| `Gravity3D` | Quasistatic | Yes | No | No |
| `GravityRefState3D` | Quasistatic | Yes | No | Yes |
| `BodyForce3D` | Quasistatic | No | Yes | No |

### 18.3 `linearelasticity/faults-2d` (5 test cases)

| Test Case | # Faults | # Materials | Notes |
|-----------|----------|-------------|-------|
| `TwoBlocksStatic` | 1 | 3 | Step slip; allowZeroResidual for TriP1 |
| `ThreeBlocksStatic` | 2 | 3 | Two parallel faults |
| `OneFaultShearNoSlip` | 1 | 2 | No slip (shear deformation); zero residual allowed |
| `TwoFaultsShearNoSlip` | 2 | 3 | No slip, two faults |
| `PlanePWave` | 1 | 2 | Dynamic fault test |

### 18.4 `incompressibleelasticity/nofaults-2d` (4 test cases)

| Test Case | Formulation | Gravity | Notes |
|-----------|-------------|---------|-------|
| `UniformPressure2D` | Quasistatic | No | Pressure + displacement |
| `UniformShear2D` | Quasistatic | No | Pure shear |
| `BodyForce2D` | Quasistatic | No | With body force |
| `Gravity2D` | Quasistatic | Yes | Gravity-driven |

### 18.5 `poroelasticity/nofaults-2d` (1 test case)

| Test Case | Formulation | # Soln Subfields | Notes |
|-----------|-------------|-----------------|-------|
| `PressureGradient` | Quasistatic | 3 or 6 | displacement, pressure, trace_strain (+ time derivatives) |

---

## 19. Summary

This design achieves:

- **Zero changes to `MMSTest`.**
- **One shared base class** (`MMSTestCase`) with all common infrastructure.
- **One mixin per equation type** (~50-80 lines each), completely self-contained.
- **One fault mixin** composable with any equation type.
- **Simple concrete test cases** that define analytical functions and delegate
  setup to mixins via one-liner hook overrides.
- **Reduced `TestCases.cc` verbosity** via a single macro per variant.
- **Open/closed principle**: new equations extend without modifying existing code.
- **No templates**, no complex metaprogramming — straightforward virtual dispatch.
- **Clear ownership model** with no double-free or use-after-free risks.
- **Extensible to multiple kinematic sources** without interface changes.
- **Works for 2D and 3D** without special handling.
