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

#include "pylith/utils/array.hh" // USES int_array, string_vector

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
pylith::meshio::Xdmf::write(const char* filenameXdfm,
			    const char* filenameHDF5)
{ // write
  std::string cellType = 0;
  int numCells = 0;
  int numCorners = 0;
  int numVertices = 0;
  int spaceDim = 0;
  int numTimeSteps = 0;
  double* timeStamps = 0;
  int numFields = 0;
  string_vector fieldName;
  string_vector fieldVectorFieldType;
  string_vector fieldType;
  int_array fieldNumPoints;
  int_array fieldFiberDim;

  if (spaceDim == 1) {
    std::cout
      << "WARNING: Xdmf grids not defined for 1-D domains.\n"
      << "Skipping creation of Xdmf file associated with HDF5 file '"
      << filenameHDF5 << "'" << std::endl;
    return;
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
    << "<!-- ============================================================ -->\n";

  _writeTimeStamps(timeStamps, numTimeSteps);

  _file
    << "<Grid Name=\"domain\" GridType=\"Uniform\">\n";
  for (int iTimeStep=0; iTimeStep < numTimeSteps; ++iTimeStep) {
    _writeGridTopology(cellType.c_str(), numCells);
    _writeGridGeometry(spaceDim);
    for (int iField=0; iField < numFields; ++iField) {
      if (2 == spaceDim && 
	  std::string("vector") == fieldVectorFieldType[iField]) {
	for (int component=0; component < spaceDim; ++component)
	  _writeGridAttributeComponent(fieldName[iField].c_str(),
				       fieldVectorFieldType[iField].c_str(),
				       fieldType[iField].c_str(),
				       numTimeSteps, fieldNumPoints[iField],
				       fieldFiberDim[iField],
				       iTimeStep, component);
      } else {
	_writeGridAttribute(fieldName[iField].c_str(),
			    fieldVectorFieldType[iField].c_str(),
			    fieldType[iField].c_str(), numTimeSteps, 
			    fieldNumPoints[iField], fieldFiberDim[iField],
			    iTimeStep);
      } // if/else
    } // for
  } // for

  _file
    << "      </Grid>\n"
    << "    </Grid>\n"
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
pylith::meshio::Xdmf::_writeTimeStamps(const double* timeStamps,
				       const int numTimeSteps)
{ // _writeTimeStamps
  if (numTimeSteps > 0) {
    assert(timeStamps);
  } // if
} // _writeTimeStamps

// ----------------------------------------------------------------------
// Write grid topology information.
void
pylith::meshio::Xdmf::_writeGridTopology(const char* cellType,
					 const int numCells)
{ // _writeGridTopology
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
pylith::meshio::Xdmf::_writeGridAttributeComponent(const char* name,
						   const char* vectorFieldType,
						   const char* center,
						   const int numTimeSteps,
						   const int numPoints,
						   const int fiberDim,
						   const int iTime,
						   const int component)
{ // _writeGridAttribute
  std::string componentName = "unknown";
  switch (component) {
  case 0:
    componentName = std::string("x-") + std::string(name);
    break;
  case 1:
    componentName = std::string("y-") + std::string(name);
    break;
  case 2:
    componentName = std::string("z-") + std::string(name);
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
    << "	   Type=\""<< vectorFieldType << "\"\n"
    << "	   Center=\"" << center << "\">\n"
    << "          <DataItem ItemType=\"HyperSlab\"\n"
    << "		    Dimensions=\"1 " << numPoints << " 1\"\n"
    << "		    Type=\"HyperSlab\">\n"
    << "            <DataItem\n"
    << "	       Dimensions=\"3 3\"\n"
    << "	       Format=\"XML\">\n"
    << "              0 0 "<< component << "\n"
    << "              1 1 1\n"
    << "              1 " << numPoints << " 1\n"
    << "	    </DataItem>\n"
    << "	    <DataItem\n"
    << "	       DataType=\"Float\" Precision=\"8\"\n"
    << "	       Dimensions=\""
    << numTimeSteps << " " << numPoints << " " << fiberDim << "\"\n"
    << "	       Format=\"HDF\">\n"
    << "	      &HeavyData;:/vertex_fields/" << name << "\n"
    << "	    </DataItem>\n"
    << "	  </DataItem>\n"
    << "	</Attribute>\n";
} // _writeGridAttribute

// ----------------------------------------------------------------------
// Write grid attribute.
void
pylith::meshio::Xdmf::_writeGridAttribute(const char* name,
					  const char* vectorFieldType,
					  const char* center,
					  const int numTimeSteps,
					  const int numPoints,
					  const int fiberDim,
					  const int iTime)
{ // _writeGridAttribute
  _file
    << "	<Attribute\n"
    << "	   Name=\"" << name << "\"\n"
    << "	   Type=\"" << vectorFieldType << "\"\n"
    << "	   Center=\"" << center << "\">\n"
    << "          <DataItem ItemType=\"HyperSlab\"\n"
    << "		    Dimensions=\"1 " << numPoints << " 1\"\n"
    << "		    Type=\"HyperSlab\">\n"
    << "            <DataItem\n"
    << "	       Dimensions=\"3 3\"\n"
    << "	       Format=\"XML\">\n"
    << "              0 0 0\n"
    << "              1 1 1\n"
    << "              1 " << numPoints << " 1\n"
    << "	    </DataItem>\n"
    << "	    <DataItem\n"
    << "	       DataType=\"Float\" Precision=\"8\"\n"
    << "	       Dimensions=\""
    << numTimeSteps << " " << numPoints << " " << fiberDim << "\"\n"
    << "	       Format=\"HDF\">\n"
    << "	      &HeavyData;:/vertex_fields/" << name << "\n"
    << "	    </DataItem>\n"
    << "	  </DataItem>\n"
    << "	</Attribute>\n";
} // _writeGridAttribute


// End of file
