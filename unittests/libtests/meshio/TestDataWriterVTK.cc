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

#include "TestDataWriterVTK.hh" // Implementation of class methods

#include "data/DataWriterVTKData.hh" // USES DataWriterVTKData

#include "pylith/meshio/DataWriterVTK.hh" // USES DataWriterVTK

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include <string.h> // USES strcmp()

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterVTK );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterVTK::setUp(void)
{ // setUp
  _data = 0;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterVTK::tearDown(void)
{ // tearDown
  _mesh = 0;
  delete _data; _data = 0;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterVTK::testConstructor(void)
{ // testConstructor
  DataWriterVTK writer;

  CPPUNIT_ASSERT(0 == writer._viewer);
  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
} // testConstructor

// ----------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestDataWriterVTK::testFilename(void)
{ // testDebug
  DataWriterVTK writer;

  const char* filename = "data.vtk";
  writer.filename(filename);
  CPPUNIT_ASSERT_EQUAL(std::string(filename), writer._filename);
} // testFilename

// ----------------------------------------------------------------------
// Test timeFormat()
void
pylith::meshio::TestDataWriterVTK::testTimeFormat(void)
{ // testTimeFormat
  DataWriterVTK writer;

  const char* format = "%4.1f";
  writer.timeFormat(format);
  CPPUNIT_ASSERT_EQUAL(std::string(format), writer._timeFormat);
} // testInterpolate

// ----------------------------------------------------------------------
// Test timeConstant()
void
pylith::meshio::TestDataWriterVTK::testTimeConstant(void)
{ // testTimeConstant
  DataWriterVTK writer;

  const double value = 4.5;
  writer.timeConstant(value);
  CPPUNIT_ASSERT_EQUAL(value, writer._timeConstant);
} // testInterpolate

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterVTK::testTimeStep(void)
{ // testTimeStep
  CPPUNIT_ASSERT(!_mesh.isNull());
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK writer;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->getDimension());

  writer.filename(_data->timestepFilename);
  writer.timeFormat(_data->timeFormat);

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  const double t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(_mesh, &cs, numTimeSteps);
    writer.openTimeStep(t, _mesh, &cs);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(_mesh, &cs, numTimeSteps, label, id);
    writer.openTimeStep(t, _mesh, &cs, label, id);
  } // else

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  writer.closeTimeStep();
  writer.close();

  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);

  checkFile(_data->timestepFilename, t, _data->timeFormat);
} // testTimeStep

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterVTK::testWriteVertexField(void)
{ // testWriteVertexField
  CPPUNIT_ASSERT(!_mesh.isNull());
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK writer;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->getDimension());

  std::vector< ALE::Obj<real_section_type> > vertexFields;
  _createVertexFields(&vertexFields);

  writer.filename(_data->vertexFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numVertexFields;

  const double t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(_mesh, &cs, numTimeSteps);
    writer.openTimeStep(t, _mesh, &cs);
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(_mesh, &cs, numTimeSteps, label, id);
    writer.openTimeStep(t, _mesh, &cs, label, id);
  } // else
  for (int i=0; i < nfields; ++i) {
    writer.writeVertexField(t, _data->vertexFieldsInfo[i].name,
			    vertexFields[i], 
			    _data->vertexFieldsInfo[i].field_type,
			    _mesh);
    CPPUNIT_ASSERT(writer._wroteVertexHeader);
    CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  } // for
  writer.closeTimeStep();
  writer.close();
  CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  
  checkFile(_data->vertexFilename, t, _data->timeFormat);
} // testWriteVertexField

// ----------------------------------------------------------------------
// Test writeCellField.
void
pylith::meshio::TestDataWriterVTK::testWriteCellField(void)
{ // testWriteCellField
  CPPUNIT_ASSERT(!_mesh.isNull());
  CPPUNIT_ASSERT(0 != _data);

  DataWriterVTK writer;

  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(_mesh->getDimension());

  std::vector< ALE::Obj<real_section_type> > cellFields;
  _createCellFields(&cellFields);

  writer.filename(_data->cellFilename);
  writer.timeFormat(_data->timeFormat);

  const int nfields = _data->numCellFields;

  const double t = _data->time;
  const int numTimeSteps = 1;
  if (0 == _data->cellsLabel) {
    writer.open(_mesh, &cs, numTimeSteps);
    writer.openTimeStep(t, _mesh, &cs);
    for (int i=0; i < nfields; ++i) {
      writer.writeCellField(t, _data->cellFieldsInfo[i].name,
                            cellFields[i], 
                            _data->cellFieldsInfo[i].field_type,
                            _mesh);
      CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } else {
    const char* label = _data->cellsLabel;
    const int id = _data->labelId;
    writer.open(_mesh, &cs, numTimeSteps, label, id);
    writer.openTimeStep(t, _mesh, &cs, label, id);
    for (int i=0; i < nfields; ++i) {
      writer.writeCellField(t, _data->cellFieldsInfo[i].name,
                            cellFields[i], 
                            _data->cellFieldsInfo[i].field_type,
                            _mesh, label, id);
      CPPUNIT_ASSERT(false == writer._wroteVertexHeader);
      CPPUNIT_ASSERT(writer._wroteCellHeader);
    } // for
  } // else
  writer.closeTimeStep();
  writer.close();
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  CPPUNIT_ASSERT(false == writer._wroteCellHeader);
  
  checkFile(_data->cellFilename, t, _data->timeFormat);
} // testWriteCellField

// ----------------------------------------------------------------------
// Test _vtkFilename.
void pylith::meshio::TestDataWriterVTK::testVtkFilename(void)
{ // testVtkFilename
  DataWriterVTK writer;

  // Append info to filename if number of time steps is 0.
  writer._numTimeSteps = 0;
  writer._filename = "output.vtk";
  CPPUNIT_ASSERT_EQUAL(std::string("output_info.vtk"), writer._vtkFilename(0.0));
		       
  // Use default normalization of 1.0, remove period from time stamp.
  writer._numTimeSteps = 100;
  writer._filename = "output.vtk";
  writer.timeFormat("%05.2f");
  CPPUNIT_ASSERT_EQUAL(std::string("output_t0230.vtk"), 
		       writer._vtkFilename(2.3));
  
  // Use normalization of 20.0, remove period from time stamp.
  writer._numTimeSteps = 100;
  writer._filename = "output.vtk";
  writer.timeFormat("%05.2f");
  writer.timeConstant(20.0);
  CPPUNIT_ASSERT_EQUAL(std::string("output_t0250.vtk"), 
		       writer._vtkFilename(50.0));
} // testVtkFilename

// ----------------------------------------------------------------------
// Create vertex fields.
void
pylith::meshio::TestDataWriterVTK::_createVertexFields(
		      std::vector< ALE::Obj<real_section_type> >* fields) const
{ // _createVertexFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(!_mesh.isNull());
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numVertexFields;

    const ALE::Obj<Mesh::label_sequence>& vertices = 
      _mesh->depthStratum(0);
    const Mesh::label_sequence::iterator verticesEnd = vertices->end();

    // Set vertex fields
    fields->resize(nfields);
    for (int i=0; i < nfields; ++i) {
      (*fields)[i] = new real_section_type(_mesh->comm(), _mesh->debug());
      (*fields)[i]->setChart(_mesh->getSieve()->getChart());
      const int fiberDim = _data->vertexFieldsInfo[i].fiber_dim;
      (*fields)[i]->setFiberDimension(vertices, fiberDim);
      _mesh->allocate((*fields)[i]);

      int ipt = 0;
      for (Mesh::label_sequence::iterator v_iter=vertices->begin();
	   v_iter != verticesEnd;
	   ++v_iter, ++ipt) {
	const double* values = &_data->vertexFields[i][ipt*fiberDim];
	(*fields)[i]->updatePoint(*v_iter, values);
      } // for
      CPPUNIT_ASSERT_EQUAL(_data->numVertices, ipt);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _createVertexFields

// ----------------------------------------------------------------------
// Create cell fields.
void
pylith::meshio::TestDataWriterVTK::_createCellFields(
		      std::vector< ALE::Obj<real_section_type> >* fields) const
{ // _createCellFields
  CPPUNIT_ASSERT(0 != fields);
  CPPUNIT_ASSERT(!_mesh.isNull());
  CPPUNIT_ASSERT(0 != _data);

  try {
    const int nfields = _data->numCellFields;

    const ALE::Obj<Mesh::label_sequence>& cells = 
      (0 == _data->cellsLabel) ? 
      _mesh->depthStratum(1) :
      _mesh->getLabelStratum(_data->cellsLabel, _data->labelId);
    const Mesh::label_sequence::iterator cellsEnd = cells->end();

    // Set cell fields
    fields->resize(nfields);
    for (int i=0; i < nfields; ++i) {
      (*fields)[i] = new real_section_type(_mesh->comm(), _mesh->debug());
      (*fields)[i]->setChart(_mesh->getSieve()->getChart());
      const int fiberDim = _data->cellFieldsInfo[i].fiber_dim;
      (*fields)[i]->setFiberDimension(cells, fiberDim);
      _mesh->allocate((*fields)[i]);

      int icell = 0;
      for (Mesh::label_sequence::iterator c_iter=cells->begin();
	   c_iter != cellsEnd;
	   ++c_iter, ++icell) {
	const double* values = &_data->cellFields[i][icell*fiberDim];
	(*fields)[i]->updatePoint(*c_iter, values);
      } // for
      CPPUNIT_ASSERT_EQUAL(_data->numCells, icell);
    } // for
  } catch (const ALE::Exception& err) {
    throw std::runtime_error(err.msg());
  } // catch
} // _createCellFields

// ----------------------------------------------------------------------
// Check VTK file against archived file.
void
pylith::meshio::TestDataWriterVTK::checkFile(const char* filenameRoot,
					     const double t,
					     const char* timeFormat)
{ // checkFile

  const std::string& fileroot = filenameRoot;

  std::ostringstream buffer;
  const int indexExt = fileroot.find(".vtk");
  // Add time stamp to filename
  char sbuffer[256];
  sprintf(sbuffer, timeFormat, t);
  std::string timestamp(sbuffer);
  const int pos = timestamp.find(".");
  if (pos != timestamp.length())
    timestamp.erase(pos, 1);
  buffer
    << std::string(fileroot, 0, indexExt) << "_t" << timestamp << ".vtk";
  
  const std::string& filename = buffer.str();
  const std::string filenameE = "data/" + filename;

  std::ifstream fileInE(filenameE.c_str());
  CPPUNIT_ASSERT(fileInE.is_open());

  std::ifstream fileIn(filename.c_str());
  CPPUNIT_ASSERT(fileIn.is_open());

  const int maxLen = 256;
  char line[maxLen];
  char lineE[maxLen];

  int i = 1;
  while(!fileInE.eof()) {
    fileInE.getline(lineE, maxLen);
    fileIn.getline(line, maxLen);
    if (0 != strcmp(line, lineE)) {
      std::cerr << "Line " << i << " of file '" << filename << "' is incorrect."
		<< std::endl;
      CPPUNIT_ASSERT(false);
    } // if
    ++i;
  } // while

  fileInE.close();
  fileIn.close();
} // checkFile


// End of file 
