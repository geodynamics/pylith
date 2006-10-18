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

void pylith::feassemble::Integrator::fillSection(const Obj<section_type>& field, const Obj<section_type>& coordinates, value_type (*f)(const value_type []))
{
  const topology_type::patch_type patch = 0;
  const section_type::chart_type& chart = field->getPatch(patch);

  for(section_type::chart_type::iterator c_iter = chart.begin(); c_iter != chart.end(); ++c_iter) {
    const section_type::value_type *coords = coordinates->restrict(patch, *c_iter);
    const section_type::value_type  value  = f(coords);

    field->updatePoint(patch, *c_iter, &value);
  }
};

void pylith::feassemble::Integrator::integrateFunction_2d(const Obj<section_type>& field, const Obj<section_type>& coordinates, value_type (*f)(const value_type []))
{
  const topology_type::patch_type               patch    = 0;
  const Obj<topology_type>&                     topology = field->getTopology();
  const Obj<topology_type::label_sequence>&     elements = topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator end      = elements->end();
  value_type detJ;

  for(topology_type::label_sequence::iterator e_iter = elements->begin(); e_iter != end; ++e_iter) {
    computeElementGeometry(coordinates, *e_iter, _v0, _Jac, PETSC_NULL, detJ);
    // Element integral
    for(int d = 0; d < _dim; d++) {
      std::cout << "v0["<<d<<"]: "<<_v0[d]<<std::endl;
    }
    for(int d = 0; d < _dim; d++) {
      for(int e = 0; e < _dim; e++) {
        std::cout << "Jac["<<d<<", "<<e<<"]: " << _Jac[d*_dim+e] << std::endl;
      }
    }
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

void pylith::feassemble::Integrator::integrateFunction_3d(const Obj<section_type>& field, const Obj<section_type>& coordinates, value_type (*f)(const value_type []))
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
      value_type xi   = _points[q*3+0] + 1.0;
      value_type eta  = _points[q*3+1] + 1.0;
      value_type zeta = _points[q*3+2] + 1.0;
      _realCoords[0] = _Jac[0]*xi + _Jac[1]*eta + _Jac[2]*zeta + _v0[0];
      _realCoords[1] = _Jac[3]*xi + _Jac[4]*eta + _Jac[5]*zeta + _v0[1];
      _realCoords[2] = _Jac[6]*xi + _Jac[7]*eta + _Jac[8]*zeta + _v0[2];
      for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
        _elemVector[i] += _basis[q*NUM_BASIS_FUNCTIONS+i]*f(_realCoords)*_weights[q]*detJ;
      }
    }
    /* Assembly */
    field->updateAdd(patch, *e_iter, _elemVector);
  }
};

void pylith::feassemble::Integrator::integrateFunction(const Obj<section_type>& field, const Obj<section_type>& coordinates, value_type (*f)(const value_type []))
{
  if (_dim == 2) {
    this->integrateFunction_2d(field, coordinates, f);
  } else if (_dim == 2) {
    this->integrateFunction_3d(field, coordinates, f);
  } else {
    throw ALE::Exception("No integration function for this dimension");
  }
};

void pylith::feassemble::Integrator::integrateLaplacianAction_2d(const Obj<section_type>& X, const Obj<section_type>& F, const Obj<section_type>& coordinates)
{
  const topology_type::patch_type               patch    = 0;
  const Obj<topology_type>&                     topology = X->getTopology();
  const Obj<topology_type::label_sequence>&     elements = topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator end      = elements->end();
  value_type detJ;

  for(topology_type::label_sequence::iterator e_iter = elements->begin(); e_iter != end; ++e_iter) {
    computeElementGeometry(coordinates, *e_iter, _v0, _Jac, _invJac, detJ);
    // Element integral
    PetscMemzero(_elemVector, NUM_BASIS_FUNCTIONS*sizeof(value_type));
    PetscMemzero(_elemMatrix, NUM_BASIS_FUNCTIONS*NUM_BASIS_FUNCTIONS*sizeof(value_type));
    for(int q = 0; q < NUM_QUADRATURE_POINTS; q++) {
      for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
        _testWork[0] = _invJac[0]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*2+0] + _invJac[2]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*2+1];
        _testWork[1] = _invJac[1]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*2+0] + _invJac[3]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*2+1];
        for(int j = 0; j < NUM_BASIS_FUNCTIONS; j++) {
          _basisWork[0] = _invJac[0]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*2+0] + _invJac[2]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*2+1];
          _basisWork[1] = _invJac[1]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*2+0] + _invJac[3]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*2+1];
          _elemMatrix[i*NUM_BASIS_FUNCTIONS+j] += (_testWork[0]*_basisWork[0] + _testWork[1]*_basisWork[1])*_weights[q]*detJ;
        }
      }
    }
    // Assembly
    const value_type *xWork = X->restrict(patch, *e_iter);
    for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
      for(int j = 0; j < NUM_BASIS_FUNCTIONS; j++) {
        _elemVector[i] += _elemMatrix[i*NUM_BASIS_FUNCTIONS+j]*xWork[j];
      }
    }
    F->updateAdd(patch, *e_iter, _elemVector);
  }
};

void pylith::feassemble::Integrator::integrateLaplacianAction_3d(const Obj<section_type>& X, const Obj<section_type>& F, const Obj<section_type>& coordinates)
{
  const topology_type::patch_type               patch    = 0;
  const Obj<topology_type>&                     topology = X->getTopology();
  const Obj<topology_type::label_sequence>&     elements = topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator end      = elements->end();
  value_type detJ;

  for(topology_type::label_sequence::iterator e_iter = elements->begin(); e_iter != end; ++e_iter) {
    computeElementGeometry(coordinates, *e_iter, _v0, _Jac, _invJac, detJ);
    // Element integral
    PetscMemzero(_elemVector, NUM_BASIS_FUNCTIONS*sizeof(value_type));
    PetscMemzero(_elemMatrix, NUM_BASIS_FUNCTIONS*NUM_BASIS_FUNCTIONS*sizeof(value_type));
    for(int q = 0; q < NUM_QUADRATURE_POINTS; q++) {
      for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
        _testWork[0] = _invJac[0]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+0] + _invJac[3]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+1] + _invJac[6]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+2];
        _testWork[1] = _invJac[1]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+0] + _invJac[4]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+1] + _invJac[7]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+2];
        _testWork[2] = _invJac[2]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+0] + _invJac[5]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+1] + _invJac[8]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+2];
        for(int j = 0; j < NUM_BASIS_FUNCTIONS; j++) {
          _basisWork[0] = _invJac[0]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+0] + _invJac[3]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+1] + _invJac[6]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+2];
          _basisWork[1] = _invJac[1]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+0] + _invJac[4]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+1] + _invJac[7]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+2];
          _basisWork[2] = _invJac[2]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+0] + _invJac[5]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+1] + _invJac[8]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+2];
          _elemMatrix[i*NUM_BASIS_FUNCTIONS+j] += (_testWork[0]*_basisWork[0] + _testWork[1]*_basisWork[1] + _testWork[2]*_basisWork[2])*_weights[q]*detJ;
        }
      }
    }
    // Assembly
    const value_type *xWork = X->restrict(patch, *e_iter);
    for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
      for(int j = 0; j < NUM_BASIS_FUNCTIONS; j++) {
        _elemVector[i] += _elemMatrix[i*NUM_BASIS_FUNCTIONS+j]*xWork[j];
      }
    }
    F->updateAdd(patch, *e_iter, _elemVector);
  }
};

void pylith::feassemble::Integrator::integrateLaplacianAction(const Obj<section_type>& X, const Obj<section_type>& F, const Obj<section_type>& coordinates)
{
  if (_dim == 2) {
    this->integrateLaplacianAction_2d(X, F, coordinates);
  } else if (_dim == 2) {
    this->integrateLaplacianAction_3d(X, F, coordinates);
  } else {
    throw ALE::Exception("No integration function for this dimension");
  }
};

void pylith::feassemble::Integrator::integrateElasticAction_3d(const Obj<section_type>& X, const Obj<section_type>& F, const Obj<section_type>& coordinates)
{
  const topology_type::patch_type               patch    = 0;
  const Obj<topology_type>&                     topology = X->getTopology();
  const Obj<topology_type::label_sequence>&     elements = topology->heightStratum(patch, 0);
  const topology_type::label_sequence::iterator end      = elements->end();
  const double                                  C        = 1.0;
  value_type detJ;

  for(topology_type::label_sequence::iterator e_iter = elements->begin(); e_iter != end; ++e_iter) {
    computeElementGeometry(coordinates, *e_iter, _v0, _Jac, _invJac, detJ);
    // Element integral
    PetscMemzero(_elemVector, NUM_BASIS_FUNCTIONS*3*sizeof(value_type));
    PetscMemzero(_elemMatrix, NUM_BASIS_FUNCTIONS*NUM_BASIS_FUNCTIONS*9*sizeof(value_type));
    for(int q = 0; q < NUM_QUADRATURE_POINTS; q++) {
      for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
        _testWork[0] = _invJac[0]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+0] + _invJac[3]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+1] + _invJac[6]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+2];
        _testWork[1] = _invJac[1]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+0] + _invJac[4]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+1] + _invJac[7]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+2];
        _testWork[2] = _invJac[2]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+0] + _invJac[5]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+1] + _invJac[8]*_basisDer[(q*NUM_BASIS_FUNCTIONS+i)*3+2];
        for(int j = 0; j < NUM_BASIS_FUNCTIONS; j++) {
          _basisWork[0] = _invJac[0]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+0] + _invJac[3]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+1] + _invJac[6]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+2];
          _basisWork[1] = _invJac[1]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+0] + _invJac[4]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+1] + _invJac[7]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+2];
          _basisWork[2] = _invJac[2]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+0] + _invJac[5]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+1] + _invJac[8]*_basisDer[(q*NUM_BASIS_FUNCTIONS+j)*3+2];
          for(int c = 0; c < 3; c++) {
            for(int d = 0; d < 3; d++) {
              for(int e = 0; e < 3; e++) {
                for(int f = 0; f < 3; f++) {
                  _elemMatrix[((i*3+c)*NUM_BASIS_FUNCTIONS + j)*3+d] += 0.25*C*(_testWork[d] + _testWork[c])*(_basisWork[f] +_basisWork[e])*_weights[q]*detJ;
                }
              }
              _elemMatrix[((i*3+c)*NUM_BASIS_FUNCTIONS + j)*3+d] += _basis[q*NUM_BASIS_FUNCTIONS+i]*_basis[q*NUM_BASIS_FUNCTIONS+j]*_weights[q]*detJ;
            }
          }
        }
      }
    }
    // Assembly
    const value_type *xWork = X->restrict(patch, *e_iter);
    for(int i = 0; i < NUM_BASIS_FUNCTIONS; i++) {
      for(int c = 0; c < 3; c++) {
        for(int j = 0; j < NUM_BASIS_FUNCTIONS; j++) {
          for(int d = 0; d < 3; d++) {
            _elemVector[i*3+c] += _elemMatrix[((i*3+c)*NUM_BASIS_FUNCTIONS+j)*3+d]*xWork[j*3+d];
          }
        }
      }
    }
    F->updateAdd(patch, *e_iter, _elemVector);
  }
};

void pylith::feassemble::Integrator::integrateElasticAction(const Obj<section_type>& X, const Obj<section_type>& F, const Obj<section_type>& coordinates)
{
  if (_dim == 2) {
    //this->integrateElasticAction_2d(X, F, coordinates);
  } else if (_dim == 2) {
    this->integrateElasticAction_3d(X, F, coordinates);
  } else {
    throw ALE::Exception("No integration function for this dimension");
  }
};
