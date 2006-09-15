#include "elemVector.hh"

void pylith::feassemble::Integrator::computeElementGeometry(const Obj<section_type>& coordinates, const point_type& e, value_type v0[], value_type J[], value_type invJ[], value_type& detJ)
{
  const topology_type::patch_type patch  = 0;
  const value_type *coords = coordinates->restrict(patch, e);
  value_type        invDet;

  for(int d = 0; d < _dim; d++) {
    v0[d] = coords[d];
  }
  for(int d = 0; d < _dim; d++) {
    for(int f = 0; f < _dim; f++) {
      J[d*_dim+f] = 0.5*(coords[(f+1)*_dim+d] - coords[0*_dim+d]);
    }
  }
  if (_dim == 1) {
    detJ = J[0];
  } else if (_dim == 2) {
    detJ = J[0]*J[3] - J[1]*J[2];
  } else if (_dim == 3) {
    detJ = J[0*3+0]*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]) +
      J[0*3+1]*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]) +
      J[0*3+2]*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]);
  }
  if (invJ) {
    invDet = 1.0/detJ;
    if (_dim == 2) {
      invJ[0] =  invDet*J[3];
      invJ[1] = -invDet*J[1];
      invJ[2] = -invDet*J[2];
      invJ[3] =  invDet*J[0];
    } else if (_dim == 3) {
      // FIX: This may be wrong
      invJ[0*3+0] = invDet*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]);
      invJ[0*3+1] = invDet*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]);
      invJ[0*3+2] = invDet*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]);
      invJ[1*3+0] = invDet*(J[0*3+1]*J[2*3+2] - J[0*3+2]*J[2*3+1]);
      invJ[1*3+1] = invDet*(J[0*3+2]*J[2*3+0] - J[0*3+0]*J[2*3+2]);
      invJ[1*3+2] = invDet*(J[0*3+0]*J[2*3+1] - J[0*3+1]*J[2*3+0]);
      invJ[2*3+0] = invDet*(J[0*3+1]*J[1*3+2] - J[0*3+2]*J[1*3+1]);
      invJ[2*3+1] = invDet*(J[0*3+2]*J[1*3+0] - J[0*3+0]*J[1*3+2]);
      invJ[2*3+2] = invDet*(J[0*3+0]*J[1*3+1] - J[0*3+1]*J[1*3+0]);
    }
  }
};

void pylith::feassemble::Integrator::integrateFunction(const Obj<section_type>& field, const Obj<section_type>& coordinates, value_type (*f)(value_type []))
{
  const topology_type::patch_type               patch    = 0;
  const Obj<topology_type>&                     topology = field->getTopology();
  const Obj<topology_type::label_sequence>&     elements = topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator end      = elements->end();
  value_type detJ;

  for(topology_type::label_sequence::iterator e_iter = elements->begin(); e_iter != end; ++e_iter) {
    computeElementGeometry(coordinates, *e_iter, _v0, _Jac, PETSC_NULL, detJ);
    // Element integral
    PetscMemzero(_elemVector, NUM_BASIS_FUNCTIONS*sizeof(value_type));
    for(int q = 0; q < NUM_QUADRATURE_POINTS; q++) {
      value_type xi  = _points[q*2+0] + 1.0;
      value_type eta = _points[q*2+1] + 1.0;
      _realCoords[0] = _Jac[0]*xi + _Jac[1]*eta + _v0[0];
      _realCoords[1] = _Jac[2]*xi + _Jac[3]*eta + _v0[1];
      for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
        _elemVector[i] += _basis[q*NUM_BASIS_FUNCTIONS+i]*f(_realCoords)*_weights[q]*detJ;
      }
    }
    /* Assembly */
    field->updateAdd(patch, *e_iter, _elemVector);
  }
};
