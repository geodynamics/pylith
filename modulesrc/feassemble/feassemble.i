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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// SWIG interface
%module feassemble

// Header files for module C++ code
%{
#include "pylith/feassemble/CellGeometry.hh"
#include "pylith/feassemble/GeometryPoint1D.hh"
#include "pylith/feassemble/GeometryPoint2D.hh"
#include "pylith/feassemble/GeometryPoint3D.hh"
#include "pylith/feassemble/GeometryLine1D.hh"
#include "pylith/feassemble/GeometryLine2D.hh"
#include "pylith/feassemble/GeometryLine3D.hh"
#include "pylith/feassemble/GeometryTri2D.hh"
#include "pylith/feassemble/GeometryTri3D.hh"
#include "pylith/feassemble/GeometryQuad2D.hh"
#include "pylith/feassemble/GeometryQuad3D.hh"
#include "pylith/feassemble/GeometryTet3D.hh"
#include "pylith/feassemble/GeometryHex3D.hh"
#include "pylith/feassemble/QuadratureRefCell.hh"

#include "pylith/topology/Mesh.hh"
#include "pylith/topology/SubMesh.hh"
#include "pylith/feassemble/Quadrature.hh"
#include "pylith/feassemble/ElasticityImplicit.hh"
#include "pylith/feassemble/ElasticityExplicit.hh"
#include "pylith/feassemble/ElasticityExplicitTri3.hh"
#include "pylith/feassemble/ElasticityExplicitTet4.hh"
#include "pylith/feassemble/ElasticityImplicitLgDeform.hh"
#include "pylith/feassemble/ElasticityExplicitLgDeform.hh"
#if defined(ENABLE_CUDA)
#include "pylith/feassemble/ElasticityImplicitCUDA.hh"
#endif

%}

%include "exception.i"
%exception {
  try {
    $action
  } catch (const ALE::Exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.message());
  } catch (const std::exception& err) {
    SWIG_exception(SWIG_RuntimeError, err.what());
  } // try/catch
 } // exception

%include "typemaps.i"
%include "../include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
import_array();
%}

%include "CellGeometry.i"
%include "GeometryPoint1D.i"
%include "GeometryPoint2D.i"
%include "GeometryPoint3D.i"
%include "GeometryLine1D.i"
%include "GeometryLine2D.i"
%include "GeometryLine3D.i"
%include "GeometryTri2D.i"
%include "GeometryTri3D.i"
%include "GeometryQuad2D.i"
%include "GeometryQuad3D.i"
%include "GeometryTet3D.i"
%include "GeometryHex3D.i"
%include "QuadratureRefCell.i"

%include "Quadrature.i"
%include "Integrator.i"
%include "IntegratorElasticity.i"
%include "ElasticityImplicit.i"
%include "ElasticityExplicit.i"
%include "ElasticityExplicitTet4.i"
%include "ElasticityExplicitTri3.i"
%include "IntegratorElasticityLgDeform.i"
%include "ElasticityImplicitLgDeform.i"
%include "ElasticityExplicitLgDeform.i"
#if defined(ENABLE_CUDA)
%include "ElasticityImplicitCUDA.i"
#endif

// Template instatiation
%template(MeshQuadrature) pylith::feassemble::Quadrature<pylith::topology::Mesh>;
%template(SubMeshQuadrature) pylith::feassemble::Quadrature<pylith::topology::SubMesh>;


// End of file

