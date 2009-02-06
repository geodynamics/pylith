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
    class GeometryPoint1D;
    class GeometryPoint2D;
    class GeometryPoint3D;
    class GeometryLine1D;
    class GeometryLine2D;
    class GeometryLine3D;
    class GeometryTri2D;
    class GeometryTri3D;
    class GeometryQuad2D;
    class GeometryQuad3D;
    class GeometryTet3D;
    class GeometryHex3D;

    class Quadrature;
    class Quadrature0D;
    class Quadrature1D;
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

  } // feassemble
} // pylith


#endif // pylith_feassemble_feassemblefwd_hh


// End of file 
