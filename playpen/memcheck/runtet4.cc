// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include <petsc.h>
#include <Python.h>

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh
#include "pylith/utils/array.hh" // USES int_array

#include "pylith/meshio/MeshIOLagrit.hh"

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D
#include "pylith/feassemble/Quadrature3D.hh" // USES Quadrature3D

#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <iostream>

int
main(int argc,
     char** argv)
{ // main
  try {
    // Initialize PETSc
    PetscErrorCode err = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    CHKERRQ(err);
    
    // Initialize Python
    Py_Initialize();

    const char* filenameGmv = "strikeslip_tet4_1000m.gmv";
    const char* filenamePset = "strikeslip_tet4_1000m.pset";

    pylith::meshio::MeshIOLagrit iohandler;
    iohandler.filenameGmv(filenameGmv);
    iohandler.filenamePset(filenamePset);
    ALE::Obj<pylith::Mesh> mesh;
    iohandler.read(&mesh);

    spatialdata::geocoords::CSCart cs;
    cs.initialize();

    pylith::feassemble::Quadrature3D quadrature;
    const int cellDim = 3;
    const int numCorners = 4;
    const int numQuadPts = 1;
    const int spaceDim = 3;
    const double verticesRef[] = {
      -1.0, -1.0, -1.0,
       1.0, -1.0, -1.0,
      -1.0,  1.0, -1.0,
      -1.0, -1.0,  1.0
    };
    const double basis[] = { 0.25, 0.25, 0.25, 0.25 };
    const double basisDeriv[] = {
      -0.5, -0.5, -0.5,
       0.5,  0.0,  0.0,
       0.0,  0.5,  0.0,
       0.0,  0.0,  0.5
    };
    const double quadPtsRef[] = { -0.5, -0.5, -0.5 };
    const double quadWts[] = { 1.33333333333  };
    quadrature.initialize(verticesRef, basis, basisDeriv, quadPtsRef, quadWts,
			  cellDim, numCorners, numQuadPts, spaceDim);

    pylith::materials::ElasticIsotropic3D materialA;
    pylith::materials::ElasticIsotropic3D materialB;
    pylith::materials::ElasticIsotropic3D materialC;
    pylith::materials::ElasticIsotropic3D materialD;

    spatialdata::spatialdb::SimpleDB db;
    spatialdata::spatialdb::SimpleIOAscii dbiohandler;
    dbiohandler.filename("mat_elastic.spatialdb");
    db.ioHandler(&dbiohandler);
    db.queryType(spatialdata::spatialdb::SimpleDB::NEAREST);
      
    materialA.db(&db);
    materialA.id(1);
    materialA.label("elastic xpos");
    materialA.initialize(mesh, &cs, &quadrature);

    materialB.db(&db);
    materialB.id(2);
    materialB.label("elastic xpos");
    materialB.initialize(mesh, &cs, &quadrature);

    materialC.db(&db);
    materialC.id(3);
    materialC.label("elastic xpos");
    materialC.initialize(mesh, &cs, &quadrature);

    materialD.db(&db);
    materialD.id(4);
    materialD.label("elastic xpos");
    materialD.initialize(mesh, &cs, &quadrature);

    // Finalize Python
    Py_Finalize();

    // Finalize PETSc
    err = PetscFinalize();
    CHKERRQ(err);
  } catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
  } catch (...) {
    abort();
  } // catch
  
  return 0;
} // main


// End of file
