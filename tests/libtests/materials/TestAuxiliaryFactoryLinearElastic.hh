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

/**
 * @file tests/libtests/materials/TestAuxiliaryFactoryLinearElastic.hh
 *
 * @brief C++ TestAuxiliaryFactoryLinearElastic object.
 *
 * C++ unit testing for AuxiliaryFactoryElastic.
 */

#if !defined(pylith_materials_testauxiliaryfactorylinearelastic_hh)
#define pylith_materials_testauxiliaryfactorylinearelastic_hh

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/materials/materialsfwd.hh" // HOLDSA AuxiliaryFactoryElastic
#include "pylith/topology/Field.hh" // HOLDSA Field::SubfieldInfo
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA Coordsys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

#include <map> // USES std::map

namespace pylith {
    namespace materials {
        class TestAuxiliaryFactoryLinearElastic;
        class TestAuxiliaryFactoryLinearElastic_Data;
    } // materials
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryLinearElastic : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestAuxiliaryFactoryLinearElastic(TestAuxiliaryFactoryLinearElastic_Data* data);

    /// Destructor.
    ~TestAuxiliaryFactoryLinearElastic(void);

    /// Test adding shear modulus, bulk modulus, and reference stress/strain subfields.
    void testAdd(void);

    /// Test setValuesFromDB().
    void testSetValuesFromDB(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /// Initialze mesh, coordinate system, auxiliary field, and factory.
    void _initialize(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    AuxiliaryFactoryElastic* _factory; ///< Test subject.
    TestAuxiliaryFactoryLinearElastic_Data* _data; ///< Test data.

    pylith::topology::Mesh* _mesh; ///< Finite-element mesh.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for test subject.

}; // class TestAuxiliaryFactoryLinearElastic

// ------------------------------------------------------------------------------------------------
class pylith::materials::TestAuxiliaryFactoryLinearElastic_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestAuxiliaryFactoryLinearElastic_Data(void);

    /// Destructor
    ~TestAuxiliaryFactoryLinearElastic_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    size_t auxDim; ///< Topological dimension of auxiliary field.

    size_t dimension; ///< Spatial dimension.
    const char* meshFilename; ///< Name of file with ASCII mesh.
    spatialdata::geocoords::CoordSys* cs; ///< Coordinate system.
    spatialdata::units::Nondimensional* normalizer; ///< Scales for nondimensionalization.

    std::map<std::string, pylith::topology::Field::SubfieldInfo> subfields;
    spatialdata::spatialdb::UserFunctionDB* auxiliaryDB; ///< Spatial database with values for solution.

}; // class TestAuxiliaryFactoryLinearElastic_Data

#endif // pylith_materials_testauxiliaryfactorylinearelastic_hh

// End of file
