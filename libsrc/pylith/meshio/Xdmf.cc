// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "Xdmf.hh" // implementation of class methods

#include "HDF5.hh" // USES HDF5

#include "pylith/utils/array.hh" // USES double_array

#include <string> // USES std::string
#include <stdexcept> // USES std::runtime_error
#include <iostream> // USES std::cout
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::meshio::Xdmf::Xdmf(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::Xdmf::~Xdmf(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Write Xdmf file associated with HDF5 file.
void
pylith::meshio::Xdmf::write(const char* filenameXdmf,
			    const char* filenameHDF5)
{ // write
  assert(filenameXdmf);
  assert(filenameHDF5);

  int numCells = 0;
  int numCorners = 0;
  std::string cellType;
  int numVertices = 0;
  int spaceDim = 0;
  double_array timeStamps;
  std::vector<FieldMetadata> fieldsMetadata;

  HDF5 h5(filenameHDF5, H5F_ACC_RDONLY);

  int ndims = 0;
  hsize_t* dims = 0;
  
  h5.getDatasetDims(&dims, &ndims, "/topology", "cells");
  assert(2 == ndims);
  numCells = dims[0];
  numCorners = dims[1];

  h5.getDatasetDims(&dims, &ndims, "/geometry", "vertices");
  assert(2 == ndims);
  numVertices = dims[0];
  spaceDim = dims[1];

  if (1 == spaceDim && 2 == numCorners)
    cellType = "Polyline";
  else if (2 == spaceDim && 3 == numCorners)
    cellType = "Triangle";
  else if (2 == spaceDim && 4 == numCorners)
    cellType = "Quadrilateral";
  else if (3 == spaceDim && 4 == numCorners)
    cellType = "Tetrahedron";
  else if (3 == spaceDim && 8 == numCorners)
    cellType = "Hexahedron";
  else {
    std::ostringstream msg;
    msg << "Unknown cell type with " << numCorners
	<< " vertices for dimension " << spaceDim << ".";
    throw std::runtime_error(msg.str());
  } // else

  // :TODO: Get time stamps

  // :TODO: Get fields metadata

  const int numFields = fieldsMetadata.size();
  for (int iField=0; iField < numFields; ++iField)
    if ("vertex_field" == fieldsMetadata[iField].domain) {
      fieldsMetadata[iField].domain = "Node";
    } else if ("cell_field" == fieldsMetadata[iField].domain) {
      fieldsMetadata[iField].domain = "Cell";
    } else {
      std::ostringstream msg;
      msg << "Unknown field type '" << fieldsMetadata[iField].domain
	  << "' for field '" 
	  << fieldsMetadata[iField].domain << "'" << std::endl;
      throw std::runtime_error(msg.str());
    } // if/else

  if (1 == spaceDim) {
    std::cout
      << "WARNING: Xdmf grids not defined for 1-D domains.\n"
      << "Skipping creation of Xdmf file associated with HDF5 file '"
      << filenameHDF5 << "'" << std::endl;
    return;
  } // if

  _file.open(filenameXdmf);
  if (!_file.is_open() || !_file.good()) {
    std::ostringstream msg;
    msg << "Could not open Xdmf file '" << filenameXdmf
	<< "' for writing metadata forHDF5 file '"
	<< filenameHDF5 << "'.\n";
    throw std::runtime_error(msg.str());
  } // if

  _file
    << "<?xml version=\"1.0\" ?>\n"
    << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
    << "<!ENTITY HeavyData \"" << filenameHDF5 << "\">\n"
    << "]>\n"
    << "\n"
    << "<Xdmf>\n"
    << "  <Domain Name=\"domain\">\n";

  _writeDomainCells(numCells, numCorners);
  _writeDomainVertices(numVertices, spaceDim);

  _file
    << "    <!-- ============================================================ -->\n";

  if (timeStamps.size() > 0) {
    _file << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" "
	  << "CollectionType=\"Temporal\">\n";

    _writeTimeStamps(timeStamps);

    const int numTimeStamps = timeStamps.size();
    for (int iTimeStep=0; iTimeStep < numTimeStamps; ++iTimeStep) {
      _file
	<< "    <Grid Name=\"domain\" GridType=\"Uniform\">\n";
      _writeGridTopology(cellType.c_str(), numCells);
      _writeGridGeometry(spaceDim);
      for (int iField=0; iField < numFields; ++iField) {
	if (2 == spaceDim && 
	    std::string("vector") == fieldsMetadata[iField].vectorFieldType) {
	  for (int component=0; component < spaceDim; ++component)
	    _writeGridAttributeComponent(fieldsMetadata[iField],
					 iTimeStep, component);
	} else {
	  _writeGridAttribute(fieldsMetadata[iField],
			      iTimeStep);
	} // if/else
      } // for
      _file << "      </Grid>\n";
    } // for
    _file << "    </Grid>\n";
  } else {
    // No time steps or fields (just the mesh).
    _file << "    <Grid Name=\"domain\" GridType=\"Uniform\">\n";
    _writeGridTopology(cellType.c_str(), numCells);
    _writeGridGeometry(spaceDim);
    _file << "    </Grid>\n";
  } // if
  
  _file
    << "  </Domain>\n"
    << "</Xdmf>\n";
  _file.close();
} // write

// ----------------------------------------------------------------------
// Write domain cell information.
void
pylith::meshio::Xdmf::_writeDomainCells(const int numCells,
					const int numCorners)
{ // _writeDomainCells
  assert(_file.is_open() && _file.good());
  _file
    << "    <DataItem Name=\"cells\"\n"
    << "	      ItemType=\"Uniform\"\n"
    << "	      Format=\"HDF\"\n" 
    << "	      NumberType=\"Float\" Precision=\"8\"\n"
    << "	      Dimensions=\"" << numCells << " " << numCorners << "\">\n"
    << "      &HeavyData;:/topology/cells\n"
    << "    </DataItem>\n";

} // _writeDomainCells

// ----------------------------------------------------------------------
// Write domain vertices information.
void
pylith::meshio::Xdmf::_writeDomainVertices(const int numVertices,
					   const int spaceDim)
{ // _writeDomainVertices
  assert(_file.is_open() && _file.good());

  _file
    << "    <DataItem Name=\"vertices\"\n"
    << "	      Format=\"HDF\"\n"
    << "	      Dimensions=\"" << numVertices << " " << spaceDim << "\">\n"
    << "      &HeavyData;:/geometry/vertices\n"
    << "    </DataItem>\n";
} // _writeDomainVertices

// ----------------------------------------------------------------------
// Write time stamps.
void
pylith::meshio::Xdmf::_writeTimeStamps(const double_array& timeStamps)
{ // _writeTimeStamps
  assert(_file.is_open() && _file.good());

  const int numTimeStamps = timeStamps.size();
  _file 
    << "      <Time TimeType=\"List\">\n"
    << "        <DataItem Format=\"XML\" NumberType=\"Float\" "
    << "Dimensions=\"" << numTimeStamps << "\">\n"
    << "        ";
  for (int i=0; i < numTimeStamps; ++i)
    _file << "  " << timeStamps[i];
  _file 
    << "\n"
    << "        </DataItem>\n"
    << "      </Time>\n";
} // _writeTimeStamps

// ----------------------------------------------------------------------
// Write grid topology information.
void
pylith::meshio::Xdmf::_writeGridTopology(const char* cellType,
					 const int numCells)
{ // _writeGridTopology
  assert(_file.is_open() && _file.good());

  _file
    << "	<Topology\n"
    << "	   TopologyType=\"Triangle\"\n"
    << "	   NumberOfElements=\"" << numCells << "\">\n"
    << "	  <DataItem Reference=\"XML\">\n"
    << "	    /Xdmf/Domain/DataItem[@Name=\"cells\"]\n"
    << "	  </DataItem>\n"
    << "	</Topology>\n";
} // _writeGridTopology

// ----------------------------------------------------------------------
// Write Grid geometry.
void
pylith::meshio::Xdmf::_writeGridGeometry(const int spaceDim)
{ // _writeGridGeometry
  assert(_file.is_open() && _file.good());
  assert(2 == spaceDim || 3 == spaceDim);

  const char* geomType = (spaceDim == 3) ? "XYZ" : "XY";

  _file
    << "	<Geometry GeometryType=\"" << geomType << "\">\n"
    << "	  <DataItem Reference=\"XML\">\n"
    << "	    /Xdmf/Domain/DataItem[@Name=\"vertices\"]\n"
    << "	  </DataItem>\n"
    << "	</Geometry>\n";
} // _writeGridGeometry

// ----------------------------------------------------------------------
// Write grid attribute.
void
pylith::meshio::Xdmf::_writeGridAttribute(const FieldMetadata& metadata,
					  const int iTime)
{ // _writeGridAttribute
  assert(_file.is_open() && _file.good());

  std::string h5FullName = "";
  if (std::string("Vertex") == metadata.domain) {
    h5FullName = std::string("/vertex_fields/") + metadata.name;
  } else if (std::string("Cell") == metadata.domain) {
    h5FullName = std::string("/cell_fields/") + metadata.name;
  } else {
    std::ostringstream msg;
    msg << "Unknown domain '" << metadata.domain << "' for field '"
	<< metadata.name << "'." << std::endl;
    throw std::runtime_error(msg.str());
  } // if/else

  _file
    << "	<Attribute\n"
    << "	   Name=\"" << metadata.name << "\"\n"
    << "	   Type=\"" << metadata.vectorFieldType << "\"\n"
    << "	   Center=\"" << metadata.domain << "\">\n"
    << "          <DataItem ItemType=\"HyperSlab\"\n"
    << "		    Dimensions=\"1 " << metadata.numPoints << " 1\"\n"
    << "		    Type=\"HyperSlab\">\n"
    << "            <DataItem\n"
    << "	       Dimensions=\"3 3\"\n"
    << "	       Format=\"XML\">\n"
    << "              0 0 0\n"
    << "              1 1 1\n"
    << "              1 " << metadata.numPoints << " 1\n"
    << "	    </DataItem>\n"
    << "	    <DataItem\n"
    << "	       DataType=\"Float\" Precision=\"8\"\n"
    << "	       Dimensions=\""
    << metadata.numTimeSteps
    << " " << metadata.numPoints
    << " " << metadata.fiberDim << "\"\n"
    << "	       Format=\"HDF\">\n"
    << "	      &HeavyData;:" << h5FullName << "\n"
    << "	    </DataItem>\n"
    << "	  </DataItem>\n"
    << "	</Attribute>\n";
} // _writeGridAttribute

// ----------------------------------------------------------------------
// Write grid attribute.
void
pylith::meshio::Xdmf::_writeGridAttributeComponent(const FieldMetadata& metadata,
						   const int iTime,
						   const int component)
{ // _writeGridAttribute
  assert(_file.is_open() && _file.good());

  std::string h5FullName = "";
  if (std::string("Vertex") == metadata.domain) {
    h5FullName = std::string("/vertex_fields/") + metadata.name;
  } else if (std::string("Cell") == metadata.domain) {
    h5FullName = std::string("/cell_fields/") + metadata.name;
  } else {
    std::ostringstream msg;
    msg << "Unknown domain '" << metadata.domain << "' for field '"
	<< metadata.name << "'." << std::endl;
    throw std::runtime_error(msg.str());
  } // if/else

  std::string componentName = "unknown";
  switch (component) {
  case 0:
    componentName = std::string("x-") + std::string(metadata.name);
    break;
  case 1:
    componentName = std::string("y-") + std::string(metadata.name);
    break;
  case 2:
    componentName = std::string("z-") + std::string(metadata.name);
    break;
  default:
    { // default
      std::ostringstream msg;
      msg << "Unknown component " << component << " while writing Xdmf file.";
      std::cerr << msg.str() << std::endl;
      assert(0);
    throw std::logic_error(msg.str());
    } // default
  } // switch

  _file
    << "	<Attribute\n"
    << "	   Name=\"" << componentName << "\"\n"
    << "	   Type=\""<< metadata.vectorFieldType << "\"\n"
    << "	   Center=\"" << metadata.domain << "\">\n"
    << "          <DataItem ItemType=\"HyperSlab\"\n"
    << "		    Dimensions=\"1 " << metadata.numPoints << " 1\"\n"
    << "		    Type=\"HyperSlab\">\n"
    << "            <DataItem\n"
    << "	       Dimensions=\"3 3\"\n"
    << "	       Format=\"XML\">\n"
    << "              0 0 "<< component << "\n"
    << "              1 1 1\n"
    << "              1 " << metadata.numPoints << " 1\n"
    << "	    </DataItem>\n"
    << "	    <DataItem\n"
    << "	       DataType=\"Float\" Precision=\"8\"\n"
    << "	       Dimensions=\""
    << metadata.numTimeSteps
    << " " << metadata.numPoints
    << " " << metadata.fiberDim << "\"\n"
    << "	       Format=\"HDF\">\n"
    << "	      &HeavyData;:" << h5FullName << "\n"
    << "	    </DataItem>\n"
    << "	  </DataItem>\n"
    << "	</Attribute>\n";
} // _writeGridAttribute


// End of file
