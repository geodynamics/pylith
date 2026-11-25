% -*- markdown -*-
% -*- coding: utf-8 -*-
%
% This code is part of PyLith, developed through the Computational Infrastructure
% for Geodynamics (https://github.com/geodynamics/pylith).
%
% Copyright (c) 2010-2026 The PyLith authors and contributors
% SPDX-License-Identifier: MIT

# MMS Test Redesign: Design Specification

## Reading Guide

This document has two parts:

- **Part I — Design (Sections 1–9).**
The architecture an implementer needs to understand before writing any code.
Read this end to end.
- **Part II — Reference & Migration (Appendices A–E).**
Interface tables, the test-case inventory, the migration plan, and per-directory build details.
Consult as needed.

### The one paragraph that explains the whole design

> **Mixins are plain helper bases with no virtual methods and no knowledge of the hierarchy.**
> Every virtual hook lives on `MMSTestCase`.
> Concrete test classes implement each hook as a one-line delegation to a mixin helper.

If that sentence is clear, the rest follows.
`MMSTestCase` defines the lifecycle and the hooks.
The mixins are dumb bags of setup logic.
The concrete class is the only place the two meet, and it meets them with one-liners.

---

# Part I — Design

## 1. Overview

This document specifies the redesign of the MMS (Method of Manufactured Solutions) test infrastructure in PyLith.
The design uses an **Equation Mixin + Shared TestCase Base** pattern to eliminate duplication across governing equations while keeping each test case simple and self-contained.

### 1.1 Goals

- Make **no changes** to the existing `MMSTest` base class (`tests/src/MMSTest.hh` / `.cc`).
- Eliminate the per-equation `Test*` classes and `*_Data` structs currently duplicated across directories.
- Keep individual test cases simple: define analytical functions, configure parameters, expose mesh/discretization variants.
- Support adding new governing equations (thermoplasticity, viscoelastic poroelasticity) without modifying existing shared infrastructure.
- Compose faults orthogonally with any equation type.
- All auxiliary fields use `spatialdata::spatialdb::UserFunctionDB` (analytical functions).
- All solution fields use PETSc kernel functions (`solution_fn` signature).

### 1.2 Constraints

- `MMSTest` remains unchanged.
- Body force and gravity remain special-case fields on the problem.
- Neumann BCs remain a separate BC class (`NeumannUserFn`).
- Single kinematic source per fault for now; the design must not preclude multiple sources.
- 3D fault directories are planned.
- Additional equations are planned: thermoplasticity and poroelasticity with viscoelastic bulk rheologies (5–7+ equation/rheology combinations total).

### 1.3 Assumptions

1. The coordinate system is always Cartesian (`spatialdata::geocoords::CSCart`).
2. Scales follow elasticity conventions via `pylith::scales::ElasticityScales`.
3. Each `tests/mmstests/<equation>/<variant>` directory uses a fixed set of mesh files in its `data/` subdirectory; test cases within a directory share meshes but vary the analytical solution.
4. Multiple cell types (Tri, Quad, Tet, Hex) and polynomial orders (P1–P4, Q1–Q4) are tested.
5. Some test cases omit certain MMS checks (e.g., skip Jacobian checks for dynamic cases).
6. The test framework is Catch2 with the existing `driver_catch2.cc` entry point.
7. The build system is GNU Autotools (`Makefile.am`).

---

## 2. Class Hierarchy

```
pylith::testing::MMSTest              (UNCHANGED — tests/src/MMSTest.hh)
    │
    └── pylith::testing::MMSTestCase  (NEW — owns lifecycle + all virtual hooks)

  Helper mixins (NOT in the hierarchy above — plain bases, no virtuals):
      pylith::testing::MaterialSetup            (NEW — shared aux-DB machinery)
          ├── pylith::testing::ElasticitySetup
          ├── pylith::testing::IncompressibleElasticitySetup
          ├── pylith::testing::PoroelasticitySetup
          └── pylith::testing::ThermoplasticitySetup   (future)
      pylith::testing::FaultKinSetup            (NEW — orthogonal to equation type)
```

`MaterialSetup` is a small shared base for the three (soon more) equation mixins.
It holds the auxiliary-field database and discretization machinery they all had duplicated.
It is still a plain helper base — no virtuals, no knowledge of `MMSTestCase`.

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

**There is no diamond.**
`MMSTestCase` inherits only from `MMSTest`.
The mixins inherit only from `MaterialSetup` (equation mixins) or nothing (`FaultKinSetup`), and none of them relate to `MMSTestCase`.
The concrete class is the sole point of multiple inheritance, and its bases are unrelated to each other.

---

## 3. How the Pieces Fit Together

The design has exactly three kinds of class.
Understanding their contract removes the need to re-read the mechanics repeatedly.

**`MMSTestCase` — the lifecycle and the hooks.**
It owns everything common to every MMS test (mesh config, coordinate system, scales, time stepping, BCs, exact-solution functions, the check bitmask).
It defines a fixed `_initialize()` sequence that calls a series of virtual hooks.
Two hooks are pure virtual (the equation *must* supply them); the rest have safe defaults (no-op or a sensible base implementation).

**Mixins — plain helper bags.**
An equation mixin (`ElasticitySetup`, etc.) knows how to build its material and rheology objects and how to register its solution subfields.
`FaultKinSetup` knows how to insert cohesive cells, configure fault fields, and register fault exact solutions.
Mixins have **no virtual methods** and know nothing about `MMSTestCase` or its hooks.
They just expose ordinary helper methods.

**Concrete test cases — the wiring.**
A concrete class inherits `MMSTestCase` plus the mixins it needs.
It defines analytical functions as private static functions, configures all parameters in its constructor, and implements each `MMSTestCase` virtual hook as a one-line delegation to a mixin helper.
Because the concrete class is the only place hooks and helpers meet, the "who calls whom" question has a single, uniform answer everywhere.

Example of the one-line delegation contract:

```cpp
void UniformStrain2D::_setupMaterials(void) override {
    _createElasticityMaterials(_spaceDim, _formulation);   // helper lives on ElasticitySetup
}
std::vector<pylith::materials::Material*>
UniformStrain2D::_getMaterialPtrs(void) override {
    return _getMaterialSetupPtrs();                        // helper lives on MaterialSetup
}
```

---

## 4. `MMSTestCase` — The Shared Base Class

### 4.1 Responsibilities

- Owns all data common to every MMS test regardless of equation type.
- Implements the fixed `_initialize()` hook sequence.
- Provides defaults for mesh I/O, nondimensionalization, problem wiring, and exact-solution registration.
- Declares the two pure-virtual hooks the equation mixin's helpers satisfy (via the concrete class), and the no-op-default hooks the fault mixin's helpers satisfy (via the concrete class).

### 4.2 Header: `tests/src/MMSTestCase.hh`

```cpp
#pragma once

#include "tests/src/MMSTest.hh"

#include "spatialdata/geocoords/CSCart.hh"
#include "pylith/scales/Scales.hh"
#include "pylith/problems/Physics.hh"
#include "pylith/topology/Field.hh"
#include "pylith/bc/bcfwd.hh"
#include "pylith/materials/materialsfwd.hh"
#include "pylith/faults/faultsfwd.hh"
#include "pylith/problems/problemsfwd.hh"
#include "spatialdata/spatialdb/spatialdbfwd.hh"

#include <memory>
#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class MMSTestCase;
    }
}

class pylith::testing::MMSTestCase : public pylith::testing::MMSTest {
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

    MMSTestCase(void);
    virtual ~MMSTestCase(void);

    /** Query whether a specific check is enabled.
     * @param[in] check Bitmask value from TestCheck enum.
     * @returns true if the check is enabled.
     */
    bool checkEnabled(unsigned int check) const;

protected:

    // --- Fixed lifecycle (override MMSTest) ---

    /// Full initialization; calls hooks in fixed order. Overrides MMSTest::_initialize().
    void _initialize(void) override;

    /// Register exact solution with PetscDS. Overrides MMSTest::_setExactSolution().
    void _setExactSolution(void) override;

    // --- Hooks: pure virtual (equation mixin helper supplies via concrete class) ---

    /// Create and configure material objects.
    virtual void _setupMaterials(void) = 0;

    /// Register solution subfields with the SolutionFactory.
    virtual void _setupSolutionSubfields(pylith::problems::SolutionFactory& factory) = 0;

    // --- Hooks: default-provided (fault mixin helper supplies via concrete class) ---

    /// Read mesh from file. Default handles MeshIOAscii and MeshIOPetsc.
    virtual void _readMesh(void);

    /// Transform mesh topology for faults. Default: no-op.
    virtual void _setupFaultTopology(void);

    /// Nondimensionalize the mesh. Default uses _cs and _scales.
    virtual void _nondimensionalize(void);

    /// Configure fault auxiliary fields and kinematic sources. Default: no-op.
    virtual void _setupFaultFields(void);

    /// Assemble the problem: materials, BCs, faults, time config. Default provided.
    virtual void _setupProblem(void);

    /// Register exact solution on domain PetscDS. Default iterates _exactSolnFns.
    virtual void _registerExactSolution(PetscDS ds, PetscDM dm);

    /// Register exact solution on fault cohesive cells. Default: no-op.
    virtual void _registerExactSolutionFault(PetscDM dm);

    // --- Pointer accessors for problem wiring (default: empty) ---

    /// Raw non-owning material pointers. Default: empty. Concrete class delegates to mixin.
    virtual std::vector<pylith::materials::Material*> _getMaterialPtrs(void);

    /// Raw non-owning fault pointers. Default: empty. Concrete class delegates to mixin.
    virtual std::vector<pylith::faults::FaultCohesive*> _getFaultPtrs(void);

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

    // --- Gravity (owned; nullptr if unused) ---
    std::unique_ptr<spatialdata::spatialdb::GravityField> _gravityField;

    // --- Boundary conditions (owned) ---
    std::vector<std::unique_ptr<pylith::bc::BoundaryCondition>> _bcs;

    // --- Solution discretizations and exact-solution functions ---
    std::vector<pylith::topology::Field::Discretization> _solnDiscretizations;
    std::vector<solution_fn> _exactSolnFns;
    std::vector<solution_fn> _exactSolnDotFns;

    // --- Test configuration ---
    unsigned int _enabledChecks;  ///< Bitmask of TestCheck values.

    // Copyable/movable: all owned resources are in unique_ptr / vector<unique_ptr>,
    // so the compiler-generated move works and copy is implicitly deleted. No hand-written
    // special members are required (see Section 6).
};
```

### 4.3 Implementation: `tests/src/MMSTestCase.cc`

Required includes:

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
#include "spatialdata/spatialdb/GravityField.hh"
#include "catch2/catch_test_macros.hpp"
#include "petscdm.h"
#include "petscds.h"
```

Constructor / destructor:

```cpp
pylith::testing::MMSTestCase::MMSTestCase(void) :
    _useAsciiMesh(true),
    _spaceDim(2),
    _t(0.0),
    _dt(0.05),
    _formulation(pylith::problems::Physics::QUASISTATIC),
    _enabledChecks(CHECK_ALL) {}

pylith::testing::MMSTestCase::~MMSTestCase(void) = default;
// unique_ptr members clean up automatically; no manual delete loops.
```

`checkEnabled`:

```cpp
bool
pylith::testing::MMSTestCase::checkEnabled(unsigned int check) const {
    return (_enabledChecks & check) != 0;
}
```

`_initialize` — the hook sequence:

```cpp
void
pylith::testing::MMSTestCase::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    REQUIRE(_mesh);

    _readMesh();               // 1. Read mesh
    _setupFaultTopology();     // 2. Insert cohesive cells (no-op if no faults)
    _nondimensionalize();      // 3. Nondimensionalize
    _setupMaterials();         // 4. Create/configure materials (pure virtual)
    _setupFaultFields();       // 5. Configure fault aux fields (no-op if no faults)
    _setupProblem();           // 6. Wire up the problem object

    // 7. Create solution field and register subfields
    REQUIRE(!_solution);
    _solution = new pylith::topology::Field(*_mesh);
    _solution->setLabel("solution");
    pylith::problems::SolutionFactory factory(*_solution, _scales);
    _setupSolutionSubfields(factory);   // pure virtual
    _problem->setSolution(_solution);

    // 8. Delegate to MMSTest base for PETSc TS/SNES setup (calls _setExactSolution()).
    pylith::testing::MMSTest::_initialize();

    PYLITH_METHOD_END;
}
```

Defaults:

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

void pylith::testing::MMSTestCase::_setupFaultTopology(void) {}

void
pylith::testing::MMSTestCase::_nondimensionalize(void) {
    _mesh->setCoordSys(&_cs);
    pylith::topology::MeshOps::nondimensionalize(_mesh, _scales);
}

void pylith::testing::MMSTestCase::_setupFaultFields(void) {}

std::vector<pylith::materials::Material*>
pylith::testing::MMSTestCase::_getMaterialPtrs(void) { return {}; }

std::vector<pylith::faults::FaultCohesive*>
pylith::testing::MMSTestCase::_getFaultPtrs(void) { return {}; }

void pylith::testing::MMSTestCase::_registerExactSolutionFault(PetscDM dm) {}
```

`_setupProblem` (default).
BC pointers are extracted from the `unique_ptr` vector as raw, non-owning pointers — the problem never takes ownership.

```cpp
void
pylith::testing::MMSTestCase::_setupProblem(void) {
    REQUIRE(_problem);

    _problem->setScales(_scales);
    _problem->setGravityField(_gravityField.get());

    auto materials = _getMaterialPtrs();   // from equation mixin via concrete class
    _problem->setMaterials(materials.data(), materials.size());

    auto faults = _getFaultPtrs();         // from fault mixin via concrete class
    if (!faults.empty()) {
        _problem->setInterfaces(faults.data(), faults.size());
    }

    std::vector<pylith::bc::BoundaryCondition*> bcPtrs;
    bcPtrs.reserve(_bcs.size());
    for (auto& bc : _bcs) { bcPtrs.push_back(bc.get()); }
    _problem->setBoundaryConditions(bcPtrs.data(), bcPtrs.size());

    _problem->setStartTime(_t);
    _problem->setEndTime(_t + _dt);
    _problem->setInitialTimeStep(_dt);
    _problem->setFormulation(_formulation);
}
```

`_setExactSolution` and default registration:

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

---

## 5. Mixins

### 5.1 `MaterialSetup` — Shared Auxiliary-Field Machinery

All equation mixins previously duplicated the same members and loop: a `UserFunctionDB` for material properties, a list of aux subfield names, a list of discretizations, and a routine to apply the DB and discretizations to each material.
That machinery is factored here.

**File:** `tests/src/MaterialSetup.hh`

```cpp
#pragma once

#include "pylith/materials/materialsfwd.hh"
#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "pylith/topology/Field.hh"
#include "pylith/problems/Physics.hh"

#include <memory>
#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class MaterialSetup;
    }
}

/// Plain helper base for equation mixins. No virtual methods.
class pylith::testing::MaterialSetup {
protected:

    MaterialSetup(void);
    ~MaterialSetup(void);   // non-virtual; owns nothing polymorphic

    /** Apply _matAuxDB and _matAuxDiscretizations to a single material.
     * Called by each equation mixin once per material it creates.
     * @param[in] material Material to configure (already allocated & typed).
     * @param[in] spaceDim Spatial dimension.
     */
    void _applyMaterialAux(pylith::materials::Material* material, int spaceDim);

    /** Raw non-owning pointers to all owned materials (as Material*).
     * The concrete class returns this from _getMaterialPtrs().
     */
    std::vector<pylith::materials::Material*> _getMaterialSetupPtrs(void);

    /** Take ownership of a freshly allocated, typed material.
     * Equation mixins call this after allocating & type-configuring each material.
     */
    void _ownMaterial(std::unique_ptr<pylith::materials::Material> material);

    // --- Shared data members ---
    spatialdata::spatialdb::UserFunctionDB _matAuxDB;
    std::vector<std::string> _matAuxSubfields;
    std::vector<pylith::topology::Field::Discretization> _matAuxDiscretizations;

private:
    std::vector<std::unique_ptr<pylith::materials::Material>> _materials;  ///< owned
};
```

**File:** `tests/src/MaterialSetup.cc`

```cpp
pylith::testing::MaterialSetup::MaterialSetup(void) {
    _matAuxDB.setDescription("material auxiliary field spatial database");
}

pylith::testing::MaterialSetup::~MaterialSetup(void) = default;

void
pylith::testing::MaterialSetup::_applyMaterialAux(
    pylith::materials::Material* material, int spaceDim) {
    material->setAuxiliaryFieldDB(&_matAuxDB);
    for (size_t j = 0; j < _matAuxSubfields.size(); ++j) {
        const auto& disc = _matAuxDiscretizations[j];
        material->setAuxiliarySubfieldDiscretization(
            _matAuxSubfields[j].c_str(),
            disc.basisOrder, disc.quadOrder, spaceDim,
            pylith::topology::FieldBase::DEFAULT_BASIS,
            disc.feSpace, disc.isBasisContinuous);
    }
}

void
pylith::testing::MaterialSetup::_ownMaterial(
    std::unique_ptr<pylith::materials::Material> material) {
    _materials.push_back(std::move(material));
}

std::vector<pylith::materials::Material*>
pylith::testing::MaterialSetup::_getMaterialSetupPtrs(void) {
    std::vector<pylith::materials::Material*> ptrs;
    ptrs.reserve(_materials.size());
    for (auto& m : _materials) { ptrs.push_back(m.get()); }
    return ptrs;
}
```

### 5.2 `ElasticitySetup`

Each equation mixin now only carries its material/rheology *type* knowledge; the aux-DB machinery is inherited from `MaterialSetup`.

**File:** `tests/src/ElasticitySetup.hh`

```cpp
#pragma once

#include "tests/src/MaterialSetup.hh"
#include "pylith/materials/Elasticity.hh"
#include "pylith/materials/IsotropicLinearElasticity.hh"
#include "pylith/problems/problemsfwd.hh"

#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class ElasticitySetup;
    }
}

class pylith::testing::ElasticitySetup : public pylith::testing::MaterialSetup {
protected:

    ElasticitySetup(void) = default;
    ~ElasticitySetup(void) = default;

    struct MaterialConfig {
        std::string name;       ///< e.g., "material-id=24"
        int labelValue;
        bool useBodyForce;
        bool useReferenceState;
    };

    /** Create Elasticity materials from _matConfigs and hand them to MaterialSetup. */
    void _createElasticityMaterials(int spaceDim,
                                    pylith::problems::Physics::FormulationEnum formulation);

    std::vector<MaterialConfig> _matConfigs;
    pylith::materials::IsotropicLinearElasticity _rheology;  ///< shared across materials
};
```

**File:** `tests/src/ElasticitySetup.cc`

```cpp
void
pylith::testing::ElasticitySetup::_createElasticityMaterials(
    int spaceDim,
    pylith::problems::Physics::FormulationEnum formulation) {
    for (const auto& cfg : _matConfigs) {
        auto material = std::make_unique<pylith::materials::Elasticity>();
        material->setFormulation(formulation);
        material->useBodyForce(cfg.useBodyForce);
        material->setIdentifier("elasticity");
        material->setName(cfg.name.c_str());
        material->setLabelValue(cfg.labelValue);
        material->setBulkRheology(&_rheology);
        _rheology.useReferenceState(cfg.useReferenceState);

        _applyMaterialAux(material.get(), spaceDim);   // shared machinery
        _ownMaterial(std::move(material));
    }
}
```

Solution subfields are registered directly by the concrete class in `_setupSolutionSubfields()`.
That is a one- or two-line call to `factory.addDisplacement(...)` and, for dynamic, `factory.addVelocity(...)`.
There is no separate mixin helper for this — the call is short enough that a helper would only add indirection.

### 5.3 `IncompressibleElasticitySetup` and `PoroelasticitySetup`

These have identical structure to `ElasticitySetup`, differing only in material type, rheology type, and the flags in `MaterialConfig`.
Both inherit `MaterialSetup` and reuse `_matAuxDB`, `_matAuxSubfields`, `_matAuxDiscretizations`, `_applyMaterialAux`, `_ownMaterial`, and `_getMaterialSetupPtrs`.
There are no per-mixin prefixed duplicates.

| Mixin | Material type | Rheology type | Extra `MaterialConfig` flags |
|-------|---------------|---------------|------------------------------|
| `IncompressibleElasticitySetup` | `IncompressibleElasticity` | `IsotropicLinearIncompElasticity` | `useBodyForce` |
| `PoroelasticitySetup` | `Poroelasticity` | `IsotropicLinearPoroelasticity` | `useBodyForce`, `useTensorPermeability` |

Their `_create*Materials()` bodies mirror `_createElasticityMaterials` exactly, substituting the type in `std::make_unique<...>()` and the equation-specific flag calls.

Solution-subfield registration (done by the concrete class):

| Equation | Quasistatic subfields | Dynamic subfields |
|----------|----------------------|-------------------|
| Elasticity | displacement | displacement, velocity |
| IncompressibleElasticity | displacement, pressure | displacement, pressure, velocity |
| Poroelasticity | displacement, pressure, trace_strain | + velocity, pressure_dot, trace_strain_dot |
| Thermoplasticity (future) | displacement, temperature | + velocity |

For fault tests the concrete class appends `factory.addLagrangeMultiplierFault(...)` after the domain subfields (see Section 5.4).

### 5.4 `FaultKinSetup`

**File:** `tests/src/FaultKinSetup.hh`

```cpp
#pragma once

#include "pylith/faults/FaultCohesiveKin.hh"
#include "pylith/faults/KinSrcStep.hh"
#include "pylith/topology/topologyfwd.hh"
#include "spatialdata/spatialdb/UserFunctionDB.hh"
#include "pylith/topology/Field.hh"
#include "tests/src/MMSTest.hh"   // for solution_fn typedef

#include <memory>
#include <vector>
#include <string>

namespace pylith {
    namespace testing {
        class FaultKinSetup;
    }
}

/// Plain helper base. No virtual methods; no knowledge of MMSTestCase.
class pylith::testing::FaultKinSetup {
protected:

    FaultKinSetup(void);
    ~FaultKinSetup(void);

    struct FaultConfig {
        int cohesiveLabelValue;            ///< e.g., 100
        std::string surfaceLabelName;      ///< e.g., "fault_xpos_faces"
        int surfaceLabelValue;             ///< default 1
        std::string buriedEdgesLabelName;  ///< empty if none
        int buriedEdgesLabelValue;         ///< default 0
    };

    /** Insert cohesive cells for all configured faults; replaces *mesh in place. */
    void _transformFaultTopology(pylith::topology::Mesh*& mesh);

    /** Set auxFieldDB and discretizations on each kinematic source / fault. */
    void _configureFaultFields(int spaceDim);

    /** Append Lagrange multiplier subfield(s) to solution. */
    void _addLagrangeSubfields(pylith::problems::SolutionFactory& factory);

    /** Raw non-owning fault pointers. Concrete class returns from _getFaultPtrs(). */
    std::vector<pylith::faults::FaultCohesive*> _getFaultKinPtrs(void);

    /** Register exact solution on fault cohesive-cell PetscDS.
     * @param[in] dm Solution DM.
     * @param[in] allExactSolnFns Domain + fault exact-solution functions.
     * @param[in] allExactSolnDotFns Domain + fault time-derivative functions.
     * @param[in] numDomainSubfields Count of domain (non-fault) subfields.
     */
    void _registerFaultExactSolution(
        PetscDM dm,
        const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnFns,
        const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnDotFns,
        size_t numDomainSubfields);

    // --- Data members ---
    std::vector<FaultConfig> _faultConfigs;

    /** Kinematic sources: outer = per fault, inner = per source within a fault.
     * Inner size is currently 1; expandable without interface change. */
    std::vector<std::vector<std::unique_ptr<pylith::faults::KinSrc>>> _kinSrcs;

    spatialdata::spatialdb::UserFunctionDB _faultAuxDB;
    std::vector<std::string> _faultAuxSubfields;
    std::vector<pylith::topology::Field::Discretization> _faultAuxDiscretizations;

    std::vector<pylith::topology::Field::Discretization> _lagrangeDiscretizations;
    size_t _numFaultSolnSubfields;

private:
    std::vector<std::unique_ptr<pylith::faults::FaultCohesiveKin>> _faultObjects;  ///< owned

    /// Persistent storage for rupture names passed to setEqRuptures (must outlive the call).
    std::vector<std::vector<std::string>> _ruptureNames;
};
```

**File:** `tests/src/FaultKinSetup.cc`

```cpp
pylith::testing::FaultKinSetup::FaultKinSetup(void) :
    _numFaultSolnSubfields(0) {
    _faultAuxDB.setDescription("fault auxiliary field spatial database");
}

pylith::testing::FaultKinSetup::~FaultKinSetup(void) = default;
// unique_ptr members clean up automatically.

void
pylith::testing::FaultKinSetup::_transformFaultTopology(
    pylith::topology::Mesh*& mesh) {
    REQUIRE(_faultConfigs.size() == _kinSrcs.size());

    _faultObjects.resize(_faultConfigs.size());
    _ruptureNames.resize(_faultConfigs.size());

    for (size_t i = 0; i < _faultConfigs.size(); ++i) {
        const auto& cfg = _faultConfigs[i];
        auto fault = std::make_unique<pylith::faults::FaultCohesiveKin>();
        fault->setCohesiveLabelValue(cfg.cohesiveLabelValue);
        fault->setSurfaceLabelName(cfg.surfaceLabelName.c_str());
        if (cfg.surfaceLabelValue != 0) {
            fault->setSurfaceLabelValue(cfg.surfaceLabelValue);
        }
        if (!cfg.buriedEdgesLabelName.empty()) {
            fault->setBuriedEdgesLabelName(cfg.buriedEdgesLabelName.c_str());
            fault->setBuriedEdgesLabelValue(cfg.buriedEdgesLabelValue);
        }

        // Build persistent rupture names and raw source pointers.
        auto& names = _ruptureNames[i];
        std::vector<const char*> nameCStrs;
        std::vector<pylith::faults::KinSrc*> srcPtrs;
        for (size_t j = 0; j < _kinSrcs[i].size(); ++j) {
            names.push_back("rupture_" + std::to_string(j));
        }
        for (size_t j = 0; j < _kinSrcs[i].size(); ++j) {
            nameCStrs.push_back(names[j].c_str());       // names storage is stable now
            srcPtrs.push_back(_kinSrcs[i][j].get());
        }
        fault->setEqRuptures(nameCStrs.data(), nameCStrs.size(),
                             srcPtrs.data(), srcPtrs.size());

        pylith::topology::Mesh* meshNew = fault->transformTopology(mesh);
        delete mesh;
        mesh = meshNew;

        _faultObjects[i] = std::move(fault);
    }
}

void
pylith::testing::FaultKinSetup::_configureFaultFields(int spaceDim) {
    for (size_t iFault = 0; iFault < _faultObjects.size(); ++iFault) {
        for (auto& src : _kinSrcs[iFault]) {
            src->auxFieldDB(&_faultAuxDB);
        }
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
    for (auto& f : _faultObjects) { ptrs.push_back(f.get()); }
    return ptrs;
}

void
pylith::testing::FaultKinSetup::_registerFaultExactSolution(
    PetscDM dm,
    const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnFns,
    const std::vector<pylith::testing::MMSTest::solution_fn>& allExactSolnDotFns,
    size_t numDomainSubfields) {
    REQUIRE(!_faultObjects.empty());

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

    for (size_t i = 0; i < numDomainSubfields; ++i) {
        PylithCallPetsc(PetscDSSetExactSolution(ds, i, allExactSolnFns[i], nullptr));
        if (i < allExactSolnDotFns.size() && allExactSolnDotFns[i]) {
            PylithCallPetsc(PetscDSSetExactSolutionTimeDerivative(
                ds, i, allExactSolnDotFns[i], nullptr));
        }
    }
    for (size_t i = 0; i < _numFaultSolnSubfields; ++i) {
        const size_t idx = numDomainSubfields + i;
        PylithCallPetsc(PetscDSSetExactSolution(ds, idx, allExactSolnFns[idx], nullptr));
    }
}
```

**Multiple sources per fault (future).**
`_kinSrcs` is already `vector<vector<unique_ptr<KinSrc>>>` and `_ruptureNames` is `vector<vector<string>>`, with distinct persistent names generated per source.
`setEqRuptures` already takes arrays, so extending the inner vectors requires no interface change.

---

## 6. Ownership Model

This is the single place ownership is discussed.
The rule is uniform: **owning containers hold `std::unique_ptr`, and anything handed to the problem is a raw, non-owning pointer obtained via `.get()`.**

| Object | Owner | Storage |
|--------|-------|---------|
| `Mesh`, `TimeDependent` problem, `Solution` field | `MMSTest` (base, unchanged) | raw, as today |
| `BoundaryCondition` objects | `MMSTestCase::_bcs` | `vector<unique_ptr<>>` |
| `GravityField` | `MMSTestCase::_gravityField` | `unique_ptr<>` |
| Materials (any equation) | `MaterialSetup::_materials` | `vector<unique_ptr<Material>>` |
| Rheology | equation mixin | stack member |
| `FaultCohesiveKin` objects | `FaultKinSetup::_faultObjects` | `vector<unique_ptr<>>` |
| `KinSrc` objects | `FaultKinSetup::_kinSrcs` | `vector<vector<unique_ptr<>>>` |
| `UserFunctionDB` objects | mixins | stack members |

Because every owned resource lives in a `unique_ptr` (or a vector of them), the concrete test classes are move-only by default with no hand-written special members.
The compiler generates a correct move and implicitly deletes copy.
Variant factories (Section 7) therefore return by value safely — the returned object is moved, transferring each `unique_ptr` exactly once.
There is no double-free risk and no move-semantics boilerplate to maintain.

The problem stores only non-owning pointers (set via `setMaterials` / `setInterfaces` / `setBoundaryConditions`) and never deletes them.
Destruction order (most-derived first, then bases in reverse declaration order) frees the mixins' owned objects before `~MMSTest` destroys the problem.
Since the problem is non-owning, this is safe.

---

## 7. Concrete Test Case Pattern

Each concrete test case:

1. Inherits `MMSTestCase` + the equation mixin (+ `FaultKinSetup` if it has faults).
2. Defines analytical functions as **private static** methods.
3. Configures all parameters in its **constructor**.
4. Implements each virtual hook as a **one-line delegation** to a mixin helper.
5. Exposes **static factory methods** per mesh/discretization variant, returning by value.

### 7.1 No-fault elasticity — `UniformStrain2D`

Header:

```cpp
#pragma once

#include "tests/src/MMSTestCase.hh"
#include "tests/src/ElasticitySetup.hh"

namespace pylith { class UniformStrain2D; }

class pylith::UniformStrain2D :
    public pylith::testing::MMSTestCase,
    public pylith::testing::ElasticitySetup {
public:
    UniformStrain2D(void);
    ~UniformStrain2D(void) = default;

    static UniformStrain2D TriP1(void);
    static UniformStrain2D TriP2(void);
    static UniformStrain2D TriP3(void);
    static UniformStrain2D QuadQ1(void);
    static UniformStrain2D QuadQ2(void);
    static UniformStrain2D QuadQ3(void);

protected:
    void _setupMaterials(void) override;
    void _setupSolutionSubfields(pylith::problems::SolutionFactory& factory) override;
    std::vector<pylith::materials::Material*> _getMaterialPtrs(void) override;

private:
    static double density(double x, double y);
    static double vs(double x, double y);
    static double vp(double x, double y);
    static double disp_x(double x, double y);
    static double disp_y(double x, double y);
    static PetscErrorCode solnkernel_disp(PetscInt spaceDim, PetscReal t,
                                          const PetscReal x[], PetscInt Nc,
                                          PetscScalar* s, void* ctx);
};
```

Implementation (key parts):

```cpp
static const double STRAIN_XX = 0.1;
static const double STRAIN_YY = 0.25;
static const double STRAIN_XY = 0.3;

double pylith::UniformStrain2D::density(double x, double y) { return 2500.0; }
double pylith::UniformStrain2D::vs(double x, double y) { return 3000.0; }
double pylith::UniformStrain2D::vp(double x, double y) { return sqrt(3.0) * vs(x, y); }
double pylith::UniformStrain2D::disp_x(double x, double y) { return STRAIN_XX*x + STRAIN_XY*y; }
double pylith::UniformStrain2D::disp_y(double x, double y) { return STRAIN_XY*x + STRAIN_YY*y; }

PetscErrorCode
pylith::UniformStrain2D::solnkernel_disp(PetscInt spaceDim, PetscReal t,
                                         const PetscReal x[], PetscInt Nc,
                                         PetscScalar* s, void* ctx) {
    s[0] = disp_x(x[0], x[1]);
    s[1] = disp_y(x[0], x[1]);
    return PETSC_SUCCESS;
}

pylith::UniformStrain2D::UniformStrain2D(void) {
    GenericComponent::setName("UniformStrain2D");

    _spaceDim = 2;
    _cs.setSpaceDim(_spaceDim);
    _formulation = pylith::problems::Physics::QUASISTATIC;
    _isJacobianLinear = true;
    _enabledChecks = CHECK_ALL;

    const double lengthScale = 8.0e+3;
    pylith::scales::ElasticityScales::setQuasistaticElasticity(&_scales, lengthScale);

    _matAuxDB.addValue("density", density, "kg/m**3");
    _matAuxDB.addValue("vp", vp, "m/s");
    _matAuxDB.addValue("vs", vs, "m/s");
    _matAuxDB.setCoordSys(_cs);
    _matAuxSubfields = {"density", "shear_modulus", "bulk_modulus"};
    _matAuxDiscretizations = {
        pylith::topology::Field::Discretization(0, 1),
        pylith::topology::Field::Discretization(0, 1),
        pylith::topology::Field::Discretization(0, 1),
    };
    _matConfigs = {{"material-id=24", 24, false, false}};
    _rheology.useReferenceState(false);

    static const PylithInt constrainedDOF[2] = {0, 1};
    auto bc = std::make_unique<pylith::bc::DirichletUserFn>();
    bc->setSubfieldName("displacement");
    bc->setLabelName("boundary");
    bc->setLabelValue(1);
    bc->setConstrainedDOF(constrainedDOF, 2);
    bc->setUserFn(solnkernel_disp);
    _bcs.push_back(std::move(bc));

    _exactSolnFns = {solnkernel_disp};
}

// Hooks — one-line delegations
void pylith::UniformStrain2D::_setupMaterials(void) {
    _createElasticityMaterials(_spaceDim, _formulation);
}
void pylith::UniformStrain2D::_setupSolutionSubfields(
    pylith::problems::SolutionFactory& factory) {
    factory.addDisplacement(_solnDiscretizations[0]);
}
std::vector<pylith::materials::Material*>
pylith::UniformStrain2D::_getMaterialPtrs(void) {
    return _getMaterialSetupPtrs();
}

// Variant factory (others follow the same shape)
pylith::UniformStrain2D
pylith::UniformStrain2D::TriP2(void) {
    UniformStrain2D test;
    test._meshFilename = "data/tri.msh";
    test._useAsciiMesh = false;
    test._solnDiscretizations = { pylith::topology::Field::Discretization(2, 2) };
    return test;   // moved out; unique_ptr members transfer cleanly
}
```

### 7.2 Fault test — `TwoBlocksStatic`

This adds `FaultKinSetup` and four extra one-line hook overrides.
The constructor additionally configures `_faultConfigs`, `_kinSrcs` (via `make_unique<KinSrcStep>()`), `_faultAuxDB`, `_faultAuxSubfields`, `_faultAuxDiscretizations`, and `_numFaultSolnSubfields`, and appends the Lagrange multiplier exact-solution function to `_exactSolnFns`.

```cpp
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
pylith::TwoBlocksStatic::_getMaterialPtrs(void) { return _getMaterialSetupPtrs(); }
std::vector<pylith::faults::FaultCohesive*>
pylith::TwoBlocksStatic::_getFaultPtrs(void) { return _getFaultKinPtrs(); }

// Fault tests override domain registration to pass DM as context (for cell-side detection):
void pylith::TwoBlocksStatic::_registerExactSolution(PetscDS ds, PetscDM dm) {
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

### 7.3 Dynamic formulation

Set `_formulation = DYNAMIC`, restrict `_enabledChecks` (typically `CHECK_DISCRETIZATION | CHECK_RESIDUAL`), populate `_exactSolnDotFns`, and add velocity in `_setupSolutionSubfields`.
No other change is needed.

---

## 8. `TestCases.cc` and the `MMS_TEST_CASE` Macro

**File:** `tests/src/MMSTestMacros.hh`

```cpp
#pragma once
#include "catch2/catch_test_macros.hpp"

/// Register all four MMS checks for a variant; each is skipped at runtime if disabled.
#define MMS_TEST_CASE(TestClass, Variant)                                     \
    TEST_CASE(#TestClass "::" #Variant "::testDiscretization",                \
              "[" #TestClass "][" #Variant "][discretization]") {             \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_DISCRETIZATION))          \
            test.testDiscretization();                                        \
    }                                                                         \
    TEST_CASE(#TestClass "::" #Variant "::testResidual",                      \
              "[" #TestClass "][" #Variant "][residual]") {                   \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_RESIDUAL))                \
            test.testResidual();                                              \
    }                                                                         \
    TEST_CASE(#TestClass "::" #Variant "::testJacobianTaylorSeries",          \
              "[" #TestClass "][" #Variant "][Jacobian Taylor series]") {     \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_JACOBIAN_TAYLOR))         \
            test.testJacobianTaylorSeries();                                  \
    }                                                                         \
    TEST_CASE(#TestClass "::" #Variant "::testJacobianFiniteDiff",            \
              "[" #TestClass "][" #Variant "][Jacobian finite difference]") { \
        auto test = pylith::TestClass::Variant();                             \
        if (test.checkEnabled(                                                \
                pylith::testing::MMSTestCase::CHECK_JACOBIAN_FD))             \
            test.testJacobianFiniteDiff();                                    \
    }
```

Usage:

```cpp
#include "tests/src/MMSTestMacros.hh"
#include "UniformStrain2D.hh"
// ...

MMS_TEST_CASE(UniformStrain2D, TriP1)
MMS_TEST_CASE(UniformStrain2D, TriP2)
MMS_TEST_CASE(UniformStrain2D, QuadQ1)
// ...
```

A test case that should skip a check sets `_enabledChecks` in its constructor.
For example, use `CHECK_DISCRETIZATION | CHECK_RESIDUAL | CHECK_JACOBIAN_TAYLOR` to skip the FD Jacobian.
The macro's runtime guard then skips it.

---

## 9. Design Decisions and Rationale

**Multiple inheritance over composition.**
Equation setup and fault setup are orthogonal concerns that compose naturally.
Multiple inheritance lets the concrete class call mixin helpers directly, with no forwarding or `this`-passing boilerplate.

**Virtual hooks over CRTP/templates.**
The project uses C-style C++ with minimal templates.
Virtual dispatch is familiar, costs nothing meaningful (hooks run once per test init), and avoids template instantiation overhead.

**`unique_ptr` ownership over manual `new`/`delete`.**
Storing all owned objects in `unique_ptr` gives correct move semantics for free, makes factories-return-by-value safe with zero special-member boilerplate, and eliminates double-free and use-after-free risk.
This is why the document needs only one short ownership section (Section 6) rather than the several the earlier draft required.

**Shared `MaterialSetup` base over per-mixin duplication.**
The aux-DB, subfield-name, and discretization machinery is identical across equation types, so it lives in one place.
Each equation mixin adds only its material/rheology *type*.
This removes the prefixed-duplicate members (`_incompMatAuxDB`, `_poroMatAuxDB`, …) that were pure collision-avoidance noise.

**Bitmask `_enabledChecks` over separate booleans.**
It is compact, extensible, and queried uniformly by the macro.

---

# Part II — Reference & Migration

## Appendix A — Initialization Sequence

No-fault case (`UniformStrain2D::TriP2().testResidual()`):

```
MMSTest::testResidual()                     [unchanged]
  └─ _initialize()  →  MMSTestCase::_initialize()
        ├─ _readMesh()               [default: MeshIOPetsc]
        ├─ _setupFaultTopology()     [no-op]
        ├─ _nondimensionalize()      [default]
        ├─ _setupMaterials()         [→ _createElasticityMaterials()]
        ├─ _setupFaultFields()       [no-op]
        ├─ _setupProblem()           [default: setMaterials/BCs/time via _getMaterialPtrs()]
        ├─ create _solution + SolutionFactory
        ├─ _setupSolutionSubfields() [factory.addDisplacement()]
        ├─ _problem->setSolution()
        └─ MMSTest::_initialize()    [PETSc TS/SNES; calls _setExactSolution()]
              └─ _setExactSolution()
                    ├─ _registerExactSolution(ds, dm)   [PetscDSSetExactSolution(..., nullptr)]
                    └─ _registerExactSolutionFault(dm)  [no-op]
  └─ DMTSCheckResidual(...)
```

The fault case adds, in the same positions: `_setupFaultTopology()` → `_transformFaultTopology()` (create `FaultCohesiveKin`, `setEqRuptures`, `transformTopology`); `_setupFaultFields()` → `_configureFaultFields()`; `_setupProblem()` also calls `setInterfaces(...)` from `_getFaultPtrs()`; `_setupSolutionSubfields()` also calls `_addLagrangeSubfields()`; and `_registerExactSolutionFault()` → `_registerFaultExactSolution()`.

## Appendix B — API Reference

### B.1 External APIs used

| Class | Method | Purpose |
|-------|--------|---------|
| `SolutionFactory` | `addDisplacement` / `addVelocity` | Displacement / velocity subfield |
| `SolutionFactory` | `addPressure` / `addPressureDot` | Pressure / its time derivative |
| `SolutionFactory` | `addTraceStrain` / `addTraceStrainDot` | Trace strain / its time derivative |
| `SolutionFactory` | `addTemperature` | Temperature subfield |
| `SolutionFactory` | `addLagrangeMultiplierFault` | Fault Lagrange multiplier |
| `DirichletUserFn` | `setUserFn` / `setUserFnDot` | Analytical BC function / derivative |
| `DirichletUserFn` | `setConstrainedDOF`, `setSubfieldName`, `setLabelName`, `setLabelValue` | BC targeting |
| `NeumannUserFn` | `setUserFn`, `setSubfieldName`, `setLabelName` | Boundary integration kernel |
| `FaultCohesiveKin` | `setCohesiveLabelValue`, `setSurfaceLabelName` | Fault identification |
| `FaultCohesiveKin` | `setEqRuptures`, `transformTopology` | Sources; insert cohesive cells |
| `FaultCohesiveKin` | `setAuxiliarySubfieldDiscretization` | Fault aux discretization |
| `KinSrc` | `setOriginTime`, `auxFieldDB` | Origin time; slip-param DB |
| `UserFunctionDB` | `addValue`, `setCoordSys` | Register analytical fn; coord sys |
| `Material` | `setAuxiliaryFieldDB`, `setAuxiliarySubfieldDiscretization` | Aux DB / discretization |
| `Material` | `setFormulation`, `setBulkRheology`, `setLabelValue`, `setName`, `setIdentifier` | Material config |

### B.2 PETSc signatures

```cpp
typedef PetscErrorCode (*PetscUserFieldFunc)(
    PetscInt dim, PetscReal t, const PetscReal x[],
    PetscInt Nc, PetscScalar u[], void* ctx);

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

### B.3 `Discretization` struct

```cpp
struct Discretization {
    int basisOrder;         // polynomial order of basis
    int quadOrder;          // quadrature order
    int dimension;          // topological dimension (-1 = default)
    int cellBasis;          // cell basis type (DEFAULT_BASIS = default)
    int feSpace;            // FE space (POLYNOMIAL_SPACE = default)
    bool isBasisContinuous; // continuity (true = default)
};
// Shorthand Discretization(basisOrder, quadOrder) uses defaults for the rest.
```

> Implementer note: confirm the exact field order and the argument order of `setAuxiliarySubfieldDiscretization(...)` against the current headers before wiring the loops in `MaterialSetup::_applyMaterialAux` and `FaultKinSetup::_configureFaultFields`.

## Appendix C — File Layout

Shared infrastructure (`tests/src/`), all NEW except where noted:

```
MMSTest.hh / .cc                        (UNCHANGED)
driver_catch2.cc                        (UNCHANGED)
MMSTestCase.hh / .cc                     (NEW)
MMSTestMacros.hh                         (NEW)
MaterialSetup.hh / .cc                   (NEW)
ElasticitySetup.hh / .cc                 (NEW)
IncompressibleElasticitySetup.hh / .cc   (NEW)
PoroelasticitySetup.hh / .cc             (NEW)
FaultKinSetup.hh / .cc                   (NEW)
```

Per-directory (e.g. `linearelasticity/nofaults-2d/`): one `.hh`/`.cc` per test case (REWRITTEN), a `TestCases.cc` (SIMPLIFIED via macro), `Makefile.am` (UPDATED), and `data/` (UNCHANGED).
Removed: `TestLinearElasticity`, `TestFaultKin`, `TestIncompressibleElasticity`, `TestLinearPoroelasticity` and their `*_Data` structs.

Updated `Makefile.am` pattern:

```makefile
include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_linearelasticity_nofaults2d

TESTS = \
    run_UniformStrain2D.sh \
    run_Gravity2D.sh
    # ... one per test case

check_SCRIPTS = $(TESTS)
check_PROGRAMS = $(MMS_DRIVER)

mmstest_linearelasticity_nofaults2d_SOURCES = \
    $(top_srcdir)/tests/src/MMSTest.cc \
    $(top_srcdir)/tests/src/MMSTestCase.cc \
    $(top_srcdir)/tests/src/MaterialSetup.cc \
    $(top_srcdir)/tests/src/ElasticitySetup.cc \
    $(top_srcdir)/tests/src/driver_catch2.cc \
    TestCases.cc \
    UniformStrain2D.cc
    # ... plus other test-case .cc; add FaultKinSetup.cc for fault dirs

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

dist_noinst_HEADERS = \
    UniformStrain2D.hh
    # ... one per test case
```

> Build note: the recipe lines under `run_%.sh:` **must be tab-indented**, not spaces (an Automake/Make requirement).

## Appendix D — Migration Plan

1. **Shared infrastructure.**
Implement `MMSTestCase`, `MaterialSetup`, `ElasticitySetup`, `FaultKinSetup`, and `MMSTestMacros.hh`; update `tests/src/Makefile.am`.
2. **Pilot: `linearelasticity/nofaults-2d`.**
Rewrite one test case (`UniformStrain2D`), convert `TestCases.cc` to the macro, remove the old `Test*` sources, update `Makefile.am`, and confirm all Catch2 tags/names are preserved and tests pass.
Then migrate the remaining cases in the directory.
3. **Remaining directories, one at a time.**
`linearelasticity/faults-2d` (needs `FaultKinSetup`); `linearelasticity/nofaults-3d` (same as 2D with `_spaceDim=3`); `incompressibleelasticity/nofaults-2d` (add `IncompressibleElasticitySetup`); `poroelasticity/nofaults-2d` (add `PoroelasticitySetup`).
4. **After each phase.**
Run the full suite and verify no regressions in coverage, test names, or CI tag filtering.

Adding a **new equation type** later: create `<Equation>Setup.hh/.cc` inheriting `MaterialSetup` (define material type, rheology type, `MaterialConfig`, and `_create<Equation>Materials`); write concrete test cases inheriting `MMSTestCase` + the new mixin; add `MMS_TEST_CASE` entries and `Makefile.am` sources.
No edits to `MMSTest`, `MMSTestCase`, `MaterialSetup`, or any other mixin or directory are required.

Adding a **new test case** to an existing equation: create the `.hh`/`.cc` (analytical fns + constructor config + one-line hook overrides + variant factories), add `#include` and `MMS_TEST_CASE` lines to `TestCases.cc`, and add sources/headers/`run_*.sh` to `Makefile.am`.

## Appendix E — Test-Case Inventory to Migrate

### E.1 `linearelasticity/nofaults-2d`

| Test Case | Formulation | Gravity | Body Force | Ref State | Notes |
|-----------|-------------|---------|------------|-----------|-------|
| `UniformStrain2D` | Quasistatic | No | No | No | Uniform strain; all 4 checks |
| `Gravity2D` | Quasistatic | Yes | No | No | Quadratic disp; P2+ only |
| `GravityRefState2D` | Quasistatic | Yes | No | Yes | Zero disp; skip FD Jacobian |
| `BodyForce2D` | Quasistatic | No | Yes | No | Quadratic disp; P2+ only |
| `RigidBodyAcc2D` | Dynamic | No | No | No | Discretization + residual only |
| `PlanePWave2D` | Dynamic | No | No | No | Discretization + residual only |

### E.2 `linearelasticity/nofaults-3d`

| Test Case | Formulation | Gravity | Body Force | Ref State |
|-----------|-------------|---------|------------|-----------|
| `UniformStrain3D` | Quasistatic | No | No | No |
| `Gravity3D` | Quasistatic | Yes | No | No |
| `GravityRefState3D` | Quasistatic | Yes | No | Yes |
| `BodyForce3D` | Quasistatic | No | Yes | No |

### E.3 `linearelasticity/faults-2d`

| Test Case | # Faults | # Materials | Notes |
|-----------|----------|-------------|-------|
| `TwoBlocksStatic` | 1 | 3 | Step slip; allowZeroResidual for TriP1 |
| `ThreeBlocksStatic` | 2 | 3 | Two parallel faults |
| `OneFaultShearNoSlip` | 1 | 2 | No slip; zero residual allowed |
| `TwoFaultsShearNoSlip` | 2 | 3 | No slip, two faults |
| `PlanePWave` | 1 | 2 | Dynamic fault test |

### E.4 `incompressibleelasticity/nofaults-2d`

| Test Case | Formulation | Gravity | Notes |
|-----------|-------------|---------|-------|
| `UniformPressure2D` | Quasistatic | No | Pressure + displacement |
| `UniformShear2D` | Quasistatic | No | Pure shear |
| `BodyForce2D` | Quasistatic | No | With body force |
| `Gravity2D` | Quasistatic | Yes | Gravity-driven |

### E.5 `poroelasticity/nofaults-2d`

| Test Case | Formulation | # Soln Subfields | Notes |
|-----------|-------------|-----------------|-------|
| `PressureGradient` | Quasistatic | 3 or 6 | displacement, pressure, trace_strain (+ derivatives) |

% End of file
