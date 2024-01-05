// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
// Robert L. Walker, Kegman, Inc.
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

/**
 * @file tests/libtests/sources/TestAuxiliaryFactoryPointForce.hh
 *
 * @brief C++ TestAuxiliaryFactoryPointForce object.
 *
 * C++ unit testing for AuxiliaryFactoryPointForce.
 */

#if !defined(pylith_sources_testauxiliaryfactorypointforce_hh)
#define pylith_sources_testauxiliaryfactorypointforce_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/sources/sourcesfwd.hh" // HOLDSA AuxiliaryFactoryPointForce
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

#include <map> // USES std::map

namespace pylith {
    namespace sources {
        class TestAuxiliaryFactoryPointForce;
        class TestAuxiliaryFactoryPointForce_Data;
    } // sources
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::sources::TestAuxiliaryFactoryPointForce : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestAuxiliaryFactoryPointForce(TestAuxiliaryFactoryPointForce_Data* data);

    /// Destructor.
    ~TestAuxiliaryFactoryPointForce(void);

    /// Test adding subfields.
    void testAdd(void);

    /// Test setValuesFromDB().
    void testSetValuesFromDB(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialze mesh, coordinate system, auxiliary field, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    AuxiliaryFactoryPointForce* _factory; ///< Test subject.
    TestAuxiliaryFactoryPointForce_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for test subject.

}; // class TestAuxiliaryFactoryPointForce

// ------------------------------------------------------------------------------------------------
class pylith::sources::TestAuxiliaryFactoryPointForce_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestAuxiliaryFactoryPointForce_Data(void);

    /// Destructor
    ~TestAuxiliaryFactoryPointForce_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    size_t dimension; ///< Spatial dimension.
    size_t auxDim; ///< Topological dimension of auxiliary field.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    std::map<std::string, pylith::topology::Field::SubfieldInfo> subfields;
    spatialdata::spatialdb::UserFunctionDB* auxiliaryDB; ///< Spatial database with values for solution.
    spatialdata::spatialdb::GravityField* gravityField; ///< Gravity field spatial database.

}; // class TestAuxiliaryFactoryPointForce_Data

#endif // pylith_sources_testAuxiliaryFactoryPointForce_hh

// End of file
