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

/** @file libsrc/feassemble/feassemblefwd.hh
 *
 * @brief Forward declarations for PyLith feassemble objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_feassemble_feassemblefwd_hh)
#define pylith_feassemble_feassemblefwd_hh

namespace pylith {
  namespace feassemble {

    class CellGeometry;
    class GeometryLine2D;
    class GeometryLine3D;
    class GeometryTri2D;
    class GeometryTri3D;
    class GeometryQuad2D;
    class GeometryQuad3D;
    class GeometryTet3D;
    class GeometryHex3D;

    class Quadrature;
    class QuadratureRefCell;
    class QuadratureEngine;
    class Quadrature1Din2D;
    class Quadrature1Din3D;
    class Quadrature2D;
    class Quadrature2Din3D;
    class Quadrature3D;

    class Constraint;
    class Integrator;

    class IntegratorElasticity;
    class ElasticityImplicit;
    class ElasticityExplicit;

    class ElasticityExplicitTet4;
    class ElasticityExplicitTri3;

    class IntegratorElasticityLgDeform;
    class ElasticityImplicitLgDeform;
    class ElasticityExplicitLgDeform;

    class ElasticityImplicitCUDA;

  } // feassemble
} // pylith


#endif // pylith_feassemble_feassemblefwd_hh


// End of file 
