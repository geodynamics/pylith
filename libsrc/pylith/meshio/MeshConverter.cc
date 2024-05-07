// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/meshio/MeshConverter.hh" // implementation of class methods

#include "pylith/meshio/MeshIO.hh" // USES MeshIO

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Convert mesh format.
void
pylith::meshio::MeshConverter::convert(pylith::meshio::MeshIO* writer,
                                       pylith::meshio::MeshIO* reader,
                                       const bool checkTopology) {
    assert(reader);
    assert(writer);

    pylith::topology::Mesh mesh;
    reader->read(&mesh, checkTopology);
    writer->write(&mesh);
}
