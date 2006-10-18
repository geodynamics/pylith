#include <Mesh.hh>

using ALE::Obj;

namespace pylith {
  namespace feassemble {
    class Integrator {
    public:
      typedef ALE::Mesh                 Mesh;
      typedef Mesh::topology_type       topology_type;
      typedef topology_type::point_type point_type;
      typedef Mesh::real_section_type   section_type;
      typedef section_type::value_type  value_type;
    protected:
      const int   NUM_QUADRATURE_POINTS;
      const int   NUM_BASIS_FUNCTIONS;
      const int   _dim;
      value_type *_points;
      value_type *_weights;
      value_type *_basis;
      value_type *_basisDer;
      value_type *_v0;
      value_type *_Jac;
      value_type *_invJac;
      value_type *_realCoords;
      value_type *_elemVector;
      value_type *_elemMatrix;
      value_type *_testWork;
      value_type *_basisWork;
    public:
      Integrator(int dim, int numQuadPoints, value_type points[], value_type weights[], int numBasisFuncs, value_type basis[], value_type basisDer[]) : NUM_QUADRATURE_POINTS(numQuadPoints), NUM_BASIS_FUNCTIONS(numBasisFuncs), _dim(dim) {
        _points     = new value_type[numQuadPoints*dim];
        _weights    = new value_type[numQuadPoints];
        _basis      = new value_type[numQuadPoints*numBasisFuncs];
        _basisDer   = new value_type[numQuadPoints*numBasisFuncs*dim];
        _v0         = new value_type[dim];
        _Jac        = new value_type[dim*dim];
        _invJac     = new value_type[dim*dim];
        _realCoords = new value_type[dim];
        _elemVector = new value_type[numBasisFuncs*dim];
        _elemMatrix = new value_type[numBasisFuncs*numBasisFuncs*dim*dim];
        _testWork   = new value_type[numBasisFuncs];
        _basisWork  = new value_type[numBasisFuncs];

        for(int q = 0; q < numQuadPoints; ++q) {
          for(int d = 0; d < dim; ++d) {
            _points[q*dim+d] = points[q*dim+d];
          }
          _weights[q] = weights[q];
          for(int b = 0; b < numBasisFuncs; ++b) {
            _basis[q*numBasisFuncs+b] = basis[q*numBasisFuncs+b];
            for(int d = 0; d < dim; ++d) {
              _basisDer[(q*numBasisFuncs+b)*dim+d] = basisDer[(q*numBasisFuncs+b)*dim+d];
            }
          }
        }
      };
      ~Integrator() {
        delete [] _points;
        delete [] _weights;
        delete [] _basis;
        delete [] _v0;
        delete [] _Jac;
        delete [] _invJac;
        delete [] _realCoords;
        delete [] _elemVector;
        delete [] _elemMatrix;
        delete [] _testWork;
        delete [] _basisWork;
      };
    protected:
      void integrateFunction_2d(const Obj<section_type>&, const Obj<section_type>&, value_type (*)(const value_type []));
      void integrateFunction_3d(const Obj<section_type>&, const Obj<section_type>&, value_type (*)(const value_type []));
      void integrateLaplacianAction_2d(const Obj<section_type>&, const Obj<section_type>&, const Obj<section_type>&);
      void integrateLaplacianAction_3d(const Obj<section_type>&, const Obj<section_type>&, const Obj<section_type>&);
      void integrateElasticAction_3d(const Obj<section_type>&, const Obj<section_type>&, const Obj<section_type>&);
    public:
      void computeElementGeometry(const Obj<section_type>&, const point_type&, value_type [], value_type [], value_type [], value_type&);
      void fillSection(const Obj<section_type>&, const Obj<section_type>&, value_type (*)(const value_type []));
      void integrateFunction(const Obj<section_type>&, const Obj<section_type>&, value_type (*)(const value_type []));
      void integrateLaplacianAction(const Obj<section_type>&, const Obj<section_type>&, const Obj<section_type>&);
      void integrateElasticAction(const Obj<section_type>&, const Obj<section_type>&, const Obj<section_type>&);
    };
  }
}
