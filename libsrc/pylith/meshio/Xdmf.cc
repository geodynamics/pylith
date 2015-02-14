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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "Xdmf.hh" // implementation of class methods

#include "HDF5.hh" // USES HDF5

#include "pylith/topology/FieldBase.hh" // USES FieldBase enums
#include "pylith/utils/array.hh" // USES PylithScalar_array
#include "pylith/utils/error.h" // USES PYLITH_METHOD_BEGIN/END

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
  PYLITH_METHOD_BEGIN;

  assert(filenameXdmf);
  assert(filenameHDF5);

  int numCells = 0;
  int numCorners = 0;
  std::string cellType;
  int numVertices = 0;
  int spaceDim = 0;
  scalar_array timeStamps;
  std::vector<FieldMetadata> fieldsMetadata;

  HDF5 h5(filenameHDF5, H5F_ACC_RDONLY);

  int ndims = 0;
  hsize_t* dims = 0;

  // Vertices
  h5.getDatasetDims(&dims, &ndims, "/geometry", "vertices");
  assert(2 == ndims);
  numVertices = dims[0];
  spaceDim = dims[1];
  delete[] dims; dims = 0;

  if (1 == spaceDim) {
    std::cout
      << "WARNING: Xdmf grids not defined for 1-D domains.\n"
      << "Skipping creation of Xdmf file associated with HDF5 file '"
      << filenameHDF5 << "'" << std::endl;
    return;
  } // if

  // Cells
  h5.getDatasetDims(&dims, &ndims, "/topology", "cells");
  assert(2 == ndims);
  numCells = dims[0];
  numCorners = dims[1];
  delete[] dims; dims = 0;
  int cellDim = 0;
  h5.readAttribute("/topology/cells", "cell_dim", (void*)&cellDim, H5T_NATIVE_INT);
  if (0 == cellDim && 1 == numCorners)
    cellType = "Polyvertex";
  else if (1 == cellDim && 2 == numCorners)
    cellType = "Polyline";
  else if (2 == cellDim && 3 == numCorners)
    cellType = "Triangle";
  else if (2 == cellDim && 4 == numCorners)
    cellType = "Quadrilateral";
  else if (3 == cellDim && 4 == numCorners)
    cellType = "Tetrahedron";
  else if (3 == cellDim && 8 == numCorners)
    cellType = "Hexahedron";
  else {
    std::ostringstream msg;
    msg << "Unknown cell type with " << numCorners
	<< " vertices and dimension " << cellDim << ".";
    throw std::runtime_error(msg.str());
  } // else

  // Time stamps
  _getTimeStamps(&timeStamps, h5);

  // Fields metadata
  _getFieldMetadata(&fieldsMetadata, h5);

  // Write Xdmf file.
  _file.open(filenameXdmf);
  if (!_file.is_open() || !_file.good()) {
    std::ostringstream msg;
    msg << "Could not open Xdmf file '" << filenameXdmf
	<< "' for writing metadata for HDF5 file '"
	<< filenameHDF5 << "'.\n";
    throw std::runtime_error(msg.str());
  } // if

  const std::string& h5string = std::string(filenameHDF5);
  const hsize_t indexSlash = h5string.find_last_of("/");
  const hsize_t indexPre = (indexSlash != h5string.npos) ? indexSlash+1 : 0;
  const std::string filenameHDF5Rel = h5string.substr(indexPre);

  _file
    << "<?xml version=\"1.0\" ?>\n"
    << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
    << "<!ENTITY HeavyData \"" << filenameHDF5Rel << "\">\n"
    << "]>\n"
    << "\n"
    << "<Xdmf>\n"
    << "  <Domain Name=\"domain\">\n";

  _writeDomainCells(numCells, numCorners);
  _writeDomainVertices(numVertices, spaceDim);

  _file
    << "    <!-- ============================================================ -->\n";

  const int numTimeStamps = timeStamps.size();
  const int numFields = fieldsMetadata.size();
  if (numTimeStamps > 1) {
    _file << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" "
	  << "CollectionType=\"Temporal\">\n";

    _writeTimeStamps(timeStamps);

    for (int iTimeStep=0; iTimeStep < numTimeStamps; ++iTimeStep) {
      _file
	<< "      <Grid Name=\"domain\" GridType=\"Uniform\">\n";
      _writeGridTopology(cellType.c_str(), numCells);
      _writeGridGeometry(spaceDim);

      for (int iField=0; iField < numFields; ++iField) {
	if (2 == spaceDim && 
	    (std::string("Vector") == fieldsMetadata[iField].vectorFieldType ||
	     std::string("Tensor6") == fieldsMetadata[iField].vectorFieldType  ||
	     std::string("Matrix") == fieldsMetadata[iField].vectorFieldType) ) {
	  const int fiberDim = fieldsMetadata[iField].fiberDim;
	  for (int component=0; component < fiberDim; ++component)
	    _writeGridAttributeComponent(fieldsMetadata[iField], iTimeStep, component, spaceDim);
	} else if (3 == spaceDim && 
		   (std::string("Tensor6") == fieldsMetadata[iField].vectorFieldType  ||
		    std::string("Matrix") == fieldsMetadata[iField].vectorFieldType) ) {
	  const int fiberDim = fieldsMetadata[iField].fiberDim;
	  for (int component=0; component < fiberDim; ++component)
	    _writeGridAttributeComponent(fieldsMetadata[iField], iTimeStep, component, spaceDim);
	} else {
	  _writeGridAttribute(fieldsMetadata[iField],
			      iTimeStep);
	} // if/else
      } // for
      _file << "      </Grid>\n";
    } // for
    _file << "    </Grid>\n";
  } else {
    // One time step or no time steps (just the mesh).
    _file << "    <Grid Name=\"domain\" GridType=\"Uniform\">\n";
    _writeGridTopology(cellType.c_str(), numCells);
    _writeGridGeometry(spaceDim);
    const int iTimeStep = 0;
    for (int iField=0; iField < numFields; ++iField) {
      if (2 == spaceDim && 
	  (std::string("Vector") == fieldsMetadata[iField].vectorFieldType ||
	   std::string("Tensor6") == fieldsMetadata[iField].vectorFieldType  ||
	   std::string("Matrix") == fieldsMetadata[iField].vectorFieldType) ) {
	const int fiberDim = fieldsMetadata[iField].fiberDim;
	for (int component=0; component < fiberDim; ++component)
	  _writeGridAttributeComponent(fieldsMetadata[iField], iTimeStep, component, spaceDim);
      } else if (3 == spaceDim && 
		 (std::string("Tensor6") == fieldsMetadata[iField].vectorFieldType  ||
		  std::string("Matrix") == fieldsMetadata[iField].vectorFieldType) ) {
	const int fiberDim = fieldsMetadata[iField].fiberDim;
	for (int component=0; component < fiberDim; ++component)
	  _writeGridAttributeComponent(fieldsMetadata[iField], iTimeStep, component, spaceDim);
      } else {
	_writeGridAttribute(fieldsMetadata[iField],
			    iTimeStep);
      } // if/else
    } // for
    _file << "    </Grid>\n";
  } // if
  
  _file
    << "  </Domain>\n"
    << "</Xdmf>\n";
  _file.close();

  PYLITH_METHOD_END;
} // write

// ----------------------------------------------------------------------
// Get timestamps from HDF5 file.
void
pylith::meshio::Xdmf::_getTimeStamps(scalar_array* timeStamps,
				     HDF5& h5)
{ // _getTimeStamps
  PYLITH_METHOD_BEGIN;

  assert(timeStamps);

  hsize_t* dims = 0;
  int ndims = 0;
  PylithScalar* t = 0;
  const hid_t scalartype = (sizeof(double) == sizeof(PylithScalar)) ? 
    H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

  if (h5.hasDataset("/time")) {
    h5.getDatasetDims(&dims, &ndims, "/", "time");
    assert(3 == ndims);
    const int numTimeStamps = dims[0];
    assert(1 == dims[1]);
    assert(1 == dims[2]);

    // Dataset is in chunks of 1.
    timeStamps->resize(numTimeStamps);
    for (int i=0; i < numTimeStamps; ++i) {
      h5.readDatasetChunk("/", "time", (char**)&t,
			  &dims, &ndims, i, scalartype);
      assert(3 == ndims);
      assert(1 == dims[0]);
      assert(1 == dims[1]);
      assert(1 == dims[2]);
      (*timeStamps)[i] = t[0];
    } // for
  } // if

  delete[] dims; dims = 0;
  delete[] t; t = 0;

  PYLITH_METHOD_END;
} // _getTimeStamps

// ----------------------------------------------------------------------
// Get field metadata from HDF5 file.
void
pylith::meshio::Xdmf::_getFieldMetadata(std::vector<FieldMetadata>* metadata,
					HDF5& h5)
{ // _getFieldMetadata
  PYLITH_METHOD_BEGIN;

  assert(metadata);

  string_vector fieldNames;
  hsize_t* dims = 0;
  int ndims = 0;

  int iOffset = 0;
  if (h5.hasGroup("/vertex_fields")) {
    const char* parent = "/vertex_fields";
    h5.getGroupDatasets(&fieldNames, parent);
    const int numFields = fieldNames.size();

    metadata->resize(iOffset+numFields);
    for (int i=0; i < numFields; ++i) {
      h5.getDatasetDims(&dims, &ndims, parent, fieldNames[i].c_str());
      std::string fullName = 
	std::string(parent) + std::string("/") + fieldNames[i];
      const std::string& vectorFieldString = h5.readAttribute(fullName.c_str(), "vector_field_type");
      switch(topology::FieldBase::parseVectorFieldString(vectorFieldString.c_str())) {
      case topology::FieldBase::SCALAR :
	(*metadata)[iOffset+i].vectorFieldType = "Scalar";
	break;
      case topology::FieldBase::VECTOR :
	(*metadata)[iOffset+i].vectorFieldType = "Vector";
	break;
      case topology::FieldBase::TENSOR :
	(*metadata)[iOffset+i].vectorFieldType = "Tensor6";
	break;
      default :
	(*metadata)[iOffset+i].vectorFieldType = "Matrix";
      } // switch

      (*metadata)[iOffset+i].name = fieldNames[i];
      (*metadata)[iOffset+i].domain = "Node";
      if (2 == ndims) {
	(*metadata)[iOffset+i].numTimeSteps = 1;
	(*metadata)[iOffset+i].numPoints = dims[0];
	(*metadata)[iOffset+i].fiberDim = dims[1];
      } else {
	assert(3 == ndims);
	(*metadata)[iOffset+i].numTimeSteps = dims[0];
	(*metadata)[iOffset+i].numPoints = dims[1];
	(*metadata)[iOffset+i].fiberDim = dims[2];
      } // if/else
      delete[] dims; dims = 0;
    } // for

    iOffset += fieldNames.size();
  } // if
  if (h5.hasGroup("/cell_fields")) {
    const char* parent = "/cell_fields";
    h5.getGroupDatasets(&fieldNames, parent);
    const int numFields = fieldNames.size();

    metadata->resize(iOffset+numFields);
    for (int i=0; i < numFields; ++i) {
      h5.getDatasetDims(&dims, &ndims, parent, fieldNames[i].c_str());
      std::string fullName =
	std::string(parent) + std::string("/") + fieldNames[i];
      const std::string& vectorFieldString = h5.readAttribute(fullName.c_str(), "vector_field_type");
      switch(topology::FieldBase::parseVectorFieldString(vectorFieldString.c_str())) {
      case topology::FieldBase::SCALAR :
	(*metadata)[iOffset+i].vectorFieldType = "Scalar";
	break;
      case topology::FieldBase::VECTOR :
	(*metadata)[iOffset+i].vectorFieldType = "Vector";
	break;
      case topology::FieldBase::TENSOR :
	(*metadata)[iOffset+i].vectorFieldType = "Tensor6";
	break;
      default :
	(*metadata)[iOffset+i].vectorFieldType = "Matrix";
      } // switch

      (*metadata)[iOffset+i].name = fieldNames[i];
      (*metadata)[iOffset+i].domain = "Cell";
      if (2 == ndims) {
	(*metadata)[iOffset+i].numTimeSteps = 1;
	(*metadata)[iOffset+i].numPoints = dims[0];
	(*metadata)[iOffset+i].fiberDim = dims[1];
      } else {
	assert(3 == ndims);
	(*metadata)[iOffset+i].numTimeSteps = dims[0];
	(*metadata)[iOffset+i].numPoints = dims[1];
	(*metadata)[iOffset+i].fiberDim = dims[2];
      } // if/else
      delete[] dims; dims = 0;
    } // for

    iOffset += fieldNames.size();
  } // if

#if 0 // debugging
  std::cout << "FIELD METADATA" << std::endl;
  for (int i=0; i < fieldNames.size(); ++i) {
    std::cout << "  Name: '" << (*metadata)[i].name << "'"
	      << ", type: " << (*metadata)[i].vectorFieldType
	      << ", domain: " << (*metadata)[i].domain
	      << ", #timesteps: " << (*metadata)[i].numTimeSteps
	      << ", numPoints: " << (*metadata)[i].numPoints
	      << ", fiberDim: " << (*metadata)[i].fiberDim
	      << std::endl;
  } // for
#endif

  PYLITH_METHOD_END;
} // _getFieldMetadata

// ----------------------------------------------------------------------
// Write domain cell information.
void
pylith::meshio::Xdmf::_writeDomainCells(const int numCells,
					const int numCorners)
{ // _writeDomainCells
  PYLITH_METHOD_BEGIN;

  assert(_file.is_open() && _file.good());
  _file
    << "    <DataItem Name=\"cells\"\n"
    << "	      ItemType=\"Uniform\"\n"
    << "	      Format=\"HDF\"\n" 
    << "	      NumberType=\"Float\" Precision=\"8\"\n"
    << "	      Dimensions=\"" << numCells << " " << numCorners << "\">\n"
    << "      &HeavyData;:/topology/cells\n"
    << "    </DataItem>\n";

  PYLITH_METHOD_END;
} // _writeDomainCells

// ----------------------------------------------------------------------
// Write domain vertices information.
void
pylith::meshio::Xdmf::_writeDomainVertices(const int numVertices,
					   const int spaceDim)
{ // _writeDomainVertices
  PYLITH_METHOD_BEGIN;

  assert(_file.is_open() && _file.good());

  _file
    << "    <DataItem Name=\"vertices\"\n"
    << "	      Format=\"HDF\"\n"
    << "	      Dimensions=\"" << numVertices << " " << spaceDim << "\">\n"
    << "      &HeavyData;:/geometry/vertices\n"
    << "    </DataItem>\n";

  PYLITH_METHOD_END;
} // _writeDomainVertices

// ----------------------------------------------------------------------
// Write time stamps.
void
pylith::meshio::Xdmf::_writeTimeStamps(const scalar_array& timeStamps)
{ // _writeTimeStamps
  PYLITH_METHOD_BEGIN;

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

  PYLITH_METHOD_END;
} // _writeTimeStamps

// ----------------------------------------------------------------------
// Write grid topology information.
void
pylith::meshio::Xdmf::_writeGridTopology(const char* cellType,
					 const int numCells)
{ // _writeGridTopology
  PYLITH_METHOD_BEGIN;

  assert(_file.is_open() && _file.good());

  _file
    << "	<Topology\n"
    << "	   TopologyType=\"" << cellType << "\"\n"
    << "	   NumberOfElements=\"" << numCells << "\">\n"
    << "	  <DataItem Reference=\"XML\">\n"
    << "	    /Xdmf/Domain/DataItem[@Name=\"cells\"]\n"
    << "	  </DataItem>\n"
    << "	</Topology>\n";

  PYLITH_METHOD_END;
} // _writeGridTopology

// ----------------------------------------------------------------------
// Write Grid geometry.
void
pylith::meshio::Xdmf::_writeGridGeometry(const int spaceDim)
{ // _writeGridGeometry
  PYLITH_METHOD_BEGIN;

  assert(_file.is_open() && _file.good());
  assert(2 == spaceDim || 3 == spaceDim);

  const char* geomType = (spaceDim == 3) ? "XYZ" : "XY";

  _file
    << "	<Geometry GeometryType=\"" << geomType << "\">\n"
    << "	  <DataItem Reference=\"XML\">\n"
    << "	    /Xdmf/Domain/DataItem[@Name=\"vertices\"]\n"
    << "	  </DataItem>\n"
    << "	</Geometry>\n";

  PYLITH_METHOD_END;
} // _writeGridGeometry

// ----------------------------------------------------------------------
// Write grid attribute.
void
pylith::meshio::Xdmf::_writeGridAttribute(const FieldMetadata& metadata,
					  const int iTime)
{ // _writeGridAttribute
  PYLITH_METHOD_BEGIN;

  assert(_file.is_open() && _file.good());

  std::string h5FullName = "";
  if (std::string("Node") == metadata.domain) {
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
    << "		    Dimensions=\"1 "
    << metadata.numPoints
    << " " << metadata.fiberDim << "\"\n"
    << "		    Type=\"HyperSlab\">\n"
    << "            <DataItem\n"
    << "	       Dimensions=\"3 3\"\n"
    << "	       Format=\"XML\">\n"
    << "              " << iTime << " 0 0\n"
    << "              1 1 1\n"
    << "              1 "
    << metadata.numPoints
    << " " << metadata.fiberDim << "\n"
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

  PYLITH_METHOD_END;
} // _writeGridAttribute

// ----------------------------------------------------------------------
// Write grid attribute.
void
pylith::meshio::Xdmf::_writeGridAttributeComponent(const FieldMetadata& metadata,
						   const int iTime,
						   const int component,
						   const int spaceDim)
{ // _writeGridAttribute
  PYLITH_METHOD_BEGIN;

  assert(_file.is_open() && _file.good());

  std::string h5FullName = "";
  if (std::string("Node") == metadata.domain) {
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
  if (std::string("Vector") == metadata.vectorFieldType) {
    switch (component) {
    case 0:
      componentName = std::string(metadata.name) + std::string("_x");
      break;
    case 1:
      componentName = std::string(metadata.name) + std::string("_y");
      break;
    case 2:
      componentName = std::string(metadata.name) + std::string("_z");
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
  } else if (std::string("Tensor6") == metadata.vectorFieldType) {
    switch (spaceDim) {
    case 2 :
      switch (component) {
      case 0:
	componentName = std::string(metadata.name) + std::string("_xx");
	break;
      case 1:
	componentName = std::string(metadata.name) + std::string("_yy");
	break;
      case 2:
	componentName = std::string(metadata.name) + std::string("_xy");
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
      break;
    case 3 :
      switch (component) {
      case 0:
	componentName = std::string(metadata.name) + std::string("_xx");
	break;
      case 1:
	componentName = std::string(metadata.name) + std::string("_yy");
	break;
      case 2:
	componentName = std::string(metadata.name) + std::string("_zz");
	break;
      case 3:
	componentName = std::string(metadata.name) + std::string("_xy");
	break;
      case 4:
	componentName = std::string(metadata.name) + std::string("_yz");
	break;
      case 5:
	componentName = std::string(metadata.name) + std::string("_xz");
	break;
      default: {
	std::ostringstream msg;
	msg << "Unknown component " << component << " while writing Xdmf file.";
	std::cerr << msg.str() << std::endl;
	assert(0);
	throw std::logic_error(msg.str());
      } // default
      } // switch
      break;
      
    default: {
	std::ostringstream msg;
	msg << "Unknown spatial dimension " << spaceDim << " while writing Xdmf file.";
	std::cerr << msg.str() << std::endl;
	assert(0);
	throw std::logic_error(msg.str());
    } // default
    } // switch
  } else if (std::string("Matrix") == metadata.vectorFieldType) {
    std::ostringstream sname;
    sname << metadata.name << "_" << component;
    componentName = sname.str();
  } else {
    std::ostringstream msg;
    msg << "Unknown vector field type " << metadata.vectorFieldType
	<< " while writing Xdmf file.";
    std::cerr << msg.str() << std::endl;
    assert(0);
    throw std::logic_error(msg.str());
  } // else
    
  _file
    << "	<Attribute\n"
    << "	   Name=\"" << componentName << "\"\n"
    << "	   Type=\"Scalar\"\n"
    << "	   Center=\"" << metadata.domain << "\">\n"
    << "          <DataItem ItemType=\"HyperSlab\"\n"
    << "		    Dimensions=\"1 "
    << metadata.numPoints
    << " 1\"\n"
    << "		    Type=\"HyperSlab\">\n"
    << "            <DataItem\n"
    << "	       Dimensions=\"3 3\"\n"
    << "	       Format=\"XML\">\n"
    << "              "
    << iTime
    << " 0 "
    << component << "\n"
    << "              1 1 1\n"
    << "              1 "
    << metadata.numPoints
    << " 1\n"
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

  PYLITH_METHOD_END;
} // _writeGridAttribute


// End of file
