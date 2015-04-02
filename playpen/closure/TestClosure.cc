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

#include "TestClosure.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields

#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

#define SEPARATE_FIELDS

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

// ----------------------------------------------------------------------
const int pylith::playpen::TestClosure::_spaceDim = 3;

// ----------------------------------------------------------------------
// Constructor
pylith::playpen::TestClosure::TestClosure(void) :
  _niterations(10)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::playpen::TestClosure::~TestClosure(void)
{ // destructor
} // destructor
  
// ----------------------------------------------------------------------
// Set number of iterations.
void
pylith::playpen::TestClosure::iterations(const long value) {
  _niterations = value;
}

// ----------------------------------------------------------------------
// Test restrictClosure().
void
pylith::playpen::TestClosure::testRestrictClosure(const pylith::topology::Mesh& mesh)
{ // testRestrictClosure

  const int spaceDim = _spaceDim;

  // Setup timing
  utils::EventLogger logger;
  logger.className("TestClosure");
  logger.initialize();
  const int stage = logger.registerStage("Test Closure");
  const int closureEvent = logger.registerEvent("closure");
  
  // Setup stuff for doing closure
  const ALE::Obj<SieveMesh>& sieveMesh = mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells =
    sieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  typedef ALE::SieveAlg<SieveMesh> SieveAlg;
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve,
      (size_t) pow(sieve->getMaxConeSize(), std::max(0, sieveMesh->depth())));
  ncV.clear();
  ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*sieve,
							       *cellsBegin, ncV);
  const int coneSize = ncV.getSize();

  // Setup coordinates visitor
  double_array coordsCell(coneSize*spaceDim);
  const ALE::Obj<RealSection>& coordsSection = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordsSection.isNull());
  topology::Mesh::RestrictVisitor coordsVisitor(*coordsSection,
						coordsCell.size(),
						&coordsCell[0]);
#if defined(SEPARATE_FIELDS)
  // Create fields
  pylith::topology::SolutionFields fields(mesh);
  fields.add("field A", "field_A");
  fields.add("field B", "field_B");
  topology::Field<topology::Mesh>& fieldA = fields.get("field A");
  fieldA.newSection(topology::FieldBase::VERTICES_FIELD, spaceDim);
  fieldA.allocate();
  fieldA.zero();
  fields.copyLayout("field A");

  // Setup field visitors
  double_array fieldACell(coneSize*spaceDim);
  const ALE::Obj<RealSection>& fieldASection = fields.get("field A").section();
  assert(!fieldASection.isNull());
  topology::Mesh::RestrictVisitor fieldAVisitor(*fieldASection,
						fieldACell.size(),
						&fieldACell[0]);
  double_array fieldBCell(coneSize*spaceDim);
  const ALE::Obj<RealSection>& fieldBSection = fields.get("field B").section();
  assert(!fieldBSection.isNull());
  topology::Mesh::RestrictVisitor fieldBVisitor(*fieldBSection,
						fieldBCell.size(),
						&fieldBCell[0]);
  
  double_array tmpCell(coneSize*spaceDim);

#else
  // Create fields
  pylith::topology::FieldsNew fields(mesh);
  fields.add("field A", "displacement", spaceDim, topology::FieldBase::VECTOR);
  fields.add("field B", "velocity", spaceDim, topology::FieldBase::VECTOR);
  fields.allocate(topology::FieldBase::VERTICES_FIELD, 2*spaceDim);

  // Create field visitors
  double_array fieldsCell(coneSize*2*spaceDim);
  const ALE::Obj<RealUniformSection>fieldsSection = fields.section();
  
  assert(!fieldsSection.isNull());
  topology::Mesh::RestrictVisitor fieldABVisitor(*fieldsSection,
						 fieldsCell.size(),
						 &fieldsCell[0]);
  
  double_array tmpCell(coneSize*2*spaceDim);
#endif
  const int dataSize = coneSize * spaceDim;

  ALE::LogStagePush(stage);
  logger.eventBegin(closureEvent);

  long count = 0;
  const long niterations = _niterations;
  for (long iter=0; iter < niterations; ++iter)
    for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	 c_iter != cellsEnd;
	 ++c_iter) {
      coordsVisitor.clear();
      sieveMesh->restrictClosure(*c_iter, coordsVisitor);
      //sieve->orientedConeOpt(*c_iter, coordsVisitor, coneSize, spaceDim);
      
#if defined(SEPARATE_FIELDS)
      fieldAVisitor.clear();
      sieveMesh->restrictClosure(*c_iter, fieldAVisitor);
      //sieve->orientedConeOpt(*c_iter, fieldAVisitor, coneSize, spaceDim);
      
      fieldBVisitor.clear();
      sieveMesh->restrictClosure(*c_iter, fieldBVisitor);
      //sieve->orientedConeOpt(*c_iter, fieldBVisitor, coneSize, spaceDim);
      
      // Perform trivial operation on fields
      //tmpCell = fieldACell + fieldBCell + coordsCell;
#else
      fieldsVisitor.clear();
      sieveMesh->restrictClosure(*c_iter, fieldsVisitor);
      
      // Perform trivial operation on fields
      //for (int i=0; i < dataSize; ++i) 
      //tmpCell[i] = fieldABCell[i] + fieldABCell[dataSize+i] + coordsCell[i];
#endif
      ++count;
      
    } // for

  logger.eventEnd(closureEvent);
  ALE::LogStagePop(stage);

  // Print stats
  PetscErrorCode ierr = 0;
  StageLog stageLog = 0;
  EventPerfLog eventLog = 0;
  ierr = PetscLogGetStageLog(&stageLog); PYLITH_CHECK_ERROR(ierr);
  ierr = StageLogGetEventPerfLog(stageLog, stage, &eventLog); PYLITH_CHECK_ERROR(ierr);
  EventPerfInfo eventInfo = eventLog->eventInfo[closureEvent];
  assert(1 == eventInfo.count);
  assert(count == cells->size() * _niterations);

#if defined(SEPARATE_FIELDS)
  const long nclosures = (1 + 2) * cells->size() * _niterations;
#else
  const long nclosures = (1 + 1) * cells->size() * _niterations;
#endif
  std::cout << "Number of cells: " << cells->size() << std::endl;
  std::cout << "Number of loops: " << _niterations << std::endl;
  std::cout << "Total time: " << eventInfo.time << std::endl;
  std::cout << "Average time per closure (" << nclosures << ") : "
	    << eventInfo.time/nclosures << std::endl;

} // testRestrictClosure


// End of file 
