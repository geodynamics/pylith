// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestMeshIOPetsc.hh" // Implementation of class methods

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace meshio {
        class TestMeshIOPetsc_ErrorCases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIOPetsc_ErrorCases {
public:

    // Data factory methods
    static TestMeshIO_Data* GmshNoEmbed2D(void);

    static TestMeshIO_Data* GmshNoEmbed3D(void);

    static TestMeshIO_Data* GmshNoSplit2D(void);

}; // TestMeshIOPetsc_ErrorCases

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestMeshIOPetsc::GmshNoEmbed2D::testReadError", "[TestMeshIOPetsc][Gmsh][Tri][binary][testReadError]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_ErrorCases::GmshNoEmbed2D()).testReadError();
}
TEST_CASE("TestMeshIOPetsc::GmshNoEmbed3D::testReadError", "[TestMeshIOPetsc][Gmsh][Tet][binary][testReadError]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_ErrorCases::GmshNoEmbed3D()).testReadError();
}
#if 0 // not a PETSc read error; adjust topology error?
TEST_CASE("TestMeshIOPetsc::GmshNoSplit2D::testReadError", "[TestMeshIOPetsc][Gmsh][Tri][binary][testReadError]") {
    pylith::meshio::TestMeshIOPetsc(pylith::meshio::TestMeshIOPetsc_ErrorCases::GmshNoSplit2D()).testReadError();
}
#endif

// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_ErrorCases::GmshNoEmbed2D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/noembed_tri.msh";

    return data;
} // GmshNoEmbed2D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_ErrorCases::GmshNoEmbed3D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/noembed_tet.msh";

    return data;
} // GmshNoEmbed2D


// ------------------------------------------------------------------------------------------------
pylith::meshio::TestMeshIO_Data*
pylith::meshio::TestMeshIOPetsc_ErrorCases::GmshNoSplit2D(void) {
    TestMeshIO_Data* data = new TestMeshIO_Data();assert(data);

    data->filename = "data/nosplit_tri.msh";

    return data;
} // GmshNoEmbed2D


// End of file
