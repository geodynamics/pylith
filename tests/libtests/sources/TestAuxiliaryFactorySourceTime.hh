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
 * @file tests/libtests/sources/TestAuxiliaryFactorySourceTime.hh
 *
 * @brief C++ TestAuxiliaryFactorySourceTime object.
 *
 * C++ unit testing for AuxiliaryFactorySourceTime.
 */

#if !defined(pylith_sources_testauxiliaryfactorysourcetime_hh)
#define pylith_sources_testauxiliaryfactorysourcetime_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/sources/sourcesfwd.hh" // HOLDSA AuxiliaryFactorySourceTime
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

#include <map> // USES std::map

namespace pylith {
    namespace sources {
        class TestAuxiliaryFactorySourceTime;
        class TestAuxiliaryFactorySourceTime_Data;
    } // sources
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::sources::TestAuxiliaryFactorySourceTime : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestAuxiliaryFactorySourceTime(TestAuxiliaryFactorySourceTime_Data* data);

    /// Destructor.
    ~TestAuxiliaryFactorySourceTime(void);

    /// Test adding density, body force, and gravity subfields.
    void testAdd(void);

    /// Test setValuesFromDB().
    void testSetValuesFromDB(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialze mesh, coordinate system, auxiliary field, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    AuxiliaryFactorySourceTime* _factory; ///< Test subject.
    TestAuxiliaryFactorySourceTime_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for test subject.

}; // class TestAuxiliaryFactorySourceTime

// ------------------------------------------------------------------------------------------------
class pylith::sources::TestAuxiliaryFactorySourceTime_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestAuxiliaryFactorySourceTime_Data(void);

    /// Destructor
    ~TestAuxiliaryFactorySourceTime_Data(void);

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

}; // class TestAuxiliaryFactorySourceTime_Data

#endif // pylith_sources_testAuxiliaryFactorySourceTime_hh

// End of file
