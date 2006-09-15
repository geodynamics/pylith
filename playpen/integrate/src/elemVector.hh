#include <Mesh.hh>

using ALE::Obj;

namespace pylith {
  namespace feassemble {
    class Integrator {
    public:
      typedef ALE::Mesh                 Mesh;
      typedef Mesh::topology_type       topology_type;
      typedef topology_type::point_type point_type;
      typedef Mesh::section_type        section_type;
      typedef section_type::value_type  value_type;
    protected:
      const int   NUM_QUADRATURE_POINTS;
      const int   NUM_BASIS_FUNCTIONS;
      const int   _dim;
      value_type *_points;
      value_type *_weights;
      value_type *_basis;
      value_type *_v0;
      value_type *_Jac;
      value_type *_invJac;
      value_type *_realCoords;
      value_type *_elemVector;
    public:
      Integrator(int dim, int numQuadPoints, value_type points[], value_type weights[], int numBasisFuncs, value_type basis[]) : NUM_QUADRATURE_POINTS(numQuadPoints), NUM_BASIS_FUNCTIONS(numBasisFuncs), _dim(dim) {
        _points     = new value_type[numQuadPoints*dim];
        _weights    = new value_type[numQuadPoints];
        _basis      = new value_type[numQuadPoints*numBasisFuncs];
        _v0         = new value_type[dim];
        _Jac        = new value_type[dim*dim];
        _invJac     = new value_type[dim*dim];
        _realCoords = new value_type[dim];
        _elemVector = new value_type[numBasisFuncs];

        for(int q = 0; q < numQuadPoints; ++q) {
          for(int d = 0; d < dim; ++d) {
            _points[q*dim+d] = points[q*dim+d];
          }
          _weights[q] = weights[q];
          for(int b = 0; b < numBasisFuncs; ++b) {
            _basis[q*numBasisFuncs+b] = basis[q*numBasisFuncs+b];
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
      };
    public:
      void computeElementGeometry(const Obj<section_type>&, const point_type&, value_type [], value_type [], value_type [], value_type&);
      void integrateFunction(const Obj<section_type>&, const Obj<section_type>&, value_type (*)(value_type []));
    };
  }
}
