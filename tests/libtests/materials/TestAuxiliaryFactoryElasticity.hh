// -*- c++ -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


#pragma once

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/materials/materialsfwd.hh" // HOLDSA AuxiliaryFactoryElasticity
#include "pylith/scales/scalesfwd.hh" // HOLDSA Scales
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB

#include <map> // USES std::map

namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryElasticity;
        class TestAuxiliaryFactoryElasticity_Data;
    } // namespace materials
} // namespace pylith

// ------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElasticity : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
  public:
    /// Constructor.
    TestAuxiliaryFactoryElasticity(TestAuxiliaryFactoryElasticity_Data* data);

    /// Destructor.
    ~TestAuxiliaryFactoryElasticity(void);

    /// Test adding density, body force, and gravity subfields.
    void testAdd(void);

    /// Test setValuesFromDB().
    void testSetValuesFromDB(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
  protected:
    /// Initialize mesh, coordinate system, auxiliary field, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
  protected:
    AuxiliaryFactoryElasticity* _factory; ///< Test subject.
    TestAuxiliaryFactoryElasticity_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for test subject.

}; // class TestAuxiliaryFactoryElasticity

// ------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryElasticity_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
  public:
    /// Constructor
    TestAuxiliaryFactoryElasticity_Data(void);

    /// Destructor
    ~TestAuxiliaryFactoryElasticity_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
  public:
    size_t dimension; ///< Spatial dimension.
    size_t auxDim; ///< Topological dimension of auxiliary field.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    pylith::scales::Scales* scales; ///< Scales for nondimensionalization.

    std::map<std::string, pylith::topology::Field::SubfieldInfo> subfields;
    spatialdata::spatialdb::UserFunctionDB* auxiliaryDB; ///< Spatial database with values for solution.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field spatial database.

}; // class TestAuxiliaryFactoryElasticity_Data


// End of file
