// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "MeshRefiner.hh" // implementation of class methods

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/meshio/DataWriter.hh" // USES DataWriter

#include <cstring> // USES strlen()
#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::topology::MeshRefiner::MeshRefiner(void)
{ // constructor
} // constructor
 
// ----------------------------------------------------------------------
// Destructor
pylith::topology::MeshRefiner::~MeshRefiner(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Write refined mesh.
void
pylith::topology::MeshRefiner::write(meshio::DataWriter* const writer,
				     const ALE::Obj<Mesh>& mesh,
				     const spatialdata::geocoords::CoordSys* cs)
{ // write
  assert(!mesh.isNull());

  const double t = 0.0;
  const int numTimeSteps = 0;
  writer->open(mesh, cs, numTimeSteps);
  writer->openTimeStep(t, mesh, cs);
  writer->closeTimeStep();
  writer->close();
} // write


// End of file 
