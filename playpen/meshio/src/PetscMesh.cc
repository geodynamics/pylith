#include "PetscMesh.hh"

ALE::PetscMesh::PetscMesh(MPI_Comm comm, 
	       int dimension, 
	       int debug) : debug(debug), dim(dimension) {
  this->setComm(comm);
  this->topology    = sieve_type(comm, debug);
  this->coordinates = field_type(comm, debug);
  this->boundary    = field_type(comm, debug);
  this->boundaries  = foliation_type(comm, debug);
  this->distributed = false;
  this->coordinates->setTopology(this->topology);
  this->boundary->setTopology(this->topology);
  this->boundaries->setTopology(this->topology);
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::getBundle"
ALE::Obj<ALE::PetscMesh::bundle_type>
ALE::PetscMesh::getBundle(const int dim) {
  ALE_LOG_EVENT_BEGIN;
  if (this->bundles.find(dim) == this->bundles.end()) {
    Obj<bundle_type> bundle = bundle_type(this->comm(), debug);
    
    // Need to globalize indices (that is what we might use the value ints for)
    std::cout << "Creating new bundle for dim " << dim << std::endl;
    bundle->setTopology(this->topology);
    bundle->setPatch(this->topology->leaves(), bundle_type::patch_type());
    bundle->setFiberDimensionByDepth(bundle_type::patch_type(), dim, 1);
    bundle->orderPatches();
    if (this->distributed) {
      bundle->createGlobalOrder();
    }
    // "element" reorder is in vertexBundle by default, and intermediate bundles could be handled by a cell tuple
    this->bundles[dim] = bundle;
  } else {
    if (this->distributed && this->bundles[dim]->getGlobalOffsets() == NULL) {
      this->bundles[dim]->createGlobalOrder();
    }
  }
  ALE_LOG_EVENT_END;
  return this->bundles[dim];
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::getField"
ALE::Obj<ALE::PetscMesh::field_type>
ALE::PetscMesh::getField(const std::string& name) {
  if (this->fields.find(name) == this->fields.end()) {
    Obj<field_type> field = field_type(this->comm(), debug);
    
    std::cout << "Creating new field " << name << std::endl;
    field->setTopology(this->topology);
    this->fields[name] = field;
  }
  return this->fields[name];
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::hasField"
bool
ALE::PetscMesh::hasField(const std::string& name) {
  return(this->fields.find(name) != this->fields.end());
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::getFields"
ALE::Obj<std::set<std::string> >
ALE::PetscMesh::getFields() {
  Obj<std::set<std::string> > names = std::set<std::string>();
  
  for(FieldContainer::iterator f_iter = this->fields.begin(); f_iter != this->fields.end(); ++f_iter) {
    names->insert(f_iter->first);
  }
  return names;
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::buildHexFaces"
// For a hex, there are 2d faces
void
ALE::PetscMesh::buildHexFaces(int dim, std::map<int, int*> *curSimplex, PointArray *boundary, point_type& simplex) {
  PointArray *faces = NULL;
  
  if (debug > 1) {std::cout << "  Building hex faces for boundary of " << simplex << " (size " << boundary->size() << "), dim " << dim << std::endl;}
  if (dim > 3) {
    throw ALE::Exception("Cannot do hexes of dimension greater than three");
  } else if (dim > 2) {
    int nodes[24] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 5, 4,
		     1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7};
    faces = new PointArray();
    
    for(int b = 0; b < 6; b++) {
      PointArray faceBoundary = PointArray();
      point_type face;
      
      for(int c = 0; c < 4; c++) {
	faceBoundary.push_back((*boundary)[nodes[b*4+c]]);
      }
      if (debug > 1) {std::cout << "    boundary point " << (*boundary)[b] << std::endl;}
      this->buildHexFaces(dim-1, curSimplex, &faceBoundary, face);
      faces->push_back(face);
    }
  } else if (dim > 1) {
    int boundarySize = (int) boundary->size();
    faces = new PointArray();
    
    for(int b = 0; b < boundarySize; b++) {
      PointArray faceBoundary = PointArray();
      point_type face;
      
      for(int c = 0; c < 2; c++) {
	faceBoundary.push_back((*boundary)[(b+c)%boundarySize]);
      }
      if (debug > 1) {std::cout << "    boundary point " << (*boundary)[b] << std::endl;}
      this->buildHexFaces(dim-1, curSimplex, &faceBoundary, face);
      faces->push_back(face);
    }
  } else {
    if (debug > 1) {std::cout << "  Just set faces to boundary in 1d" << std::endl;}
    faces = boundary;
  }
  if (debug > 1) {
    for(PointArray::iterator f_itor = faces->begin(); f_itor != faces->end(); ++f_itor) {
      std::cout << "  face point " << *f_itor << std::endl;
    }
  }
  // We always create the toplevel, so we could short circuit somehow
  // Should not have to loop here since the meet of just 2 boundary elements is an element
  PointArray::iterator f_itor = faces->begin();
  point_type           start = *f_itor;
  f_itor++;
  point_type           next = *f_itor;
  Obj<sieve_type::supportSet> preElement = this->topology->nJoin(start, next, 1);
  
  if (preElement->size() > 0) {
    simplex = *preElement->begin();
    if (debug > 1) {std::cout << "  Found old simplex " << simplex << std::endl;}
  } else {
    int color = 0;
    
    simplex = point_type(0, (*(*curSimplex)[dim])++);
    for(PointArray::iterator f_itor = faces->begin(); f_itor != faces->end(); ++f_itor) {
      this->topology->addArrow(*f_itor, simplex, color++);
    }
    if (debug > 1) {std::cout << "  Added simplex " << simplex << " dim " << dim << std::endl;}
  }
  if (dim > 1) {
    delete faces;
  }
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::buildFaces"
void
ALE::PetscMesh::buildFaces(int dim, std::map<int, int*> *curSimplex, PointArray *boundary, point_type& simplex) {
  PointArray *faces = NULL;
  
  if (debug > 1) {std::cout << "  Building faces for boundary of " << simplex << " (size " << boundary->size() << "), dim " << dim << std::endl;}
  if (dim > 1) {
    faces = new PointArray();
    
    // Use the cone construction
    for(PointArray::iterator b_itor = boundary->begin(); b_itor != boundary->end(); ++b_itor) {
      PointArray faceBoundary = PointArray();
      point_type face;
      
      for(PointArray::iterator i_itor = boundary->begin(); i_itor != boundary->end(); ++i_itor) {
	if (i_itor != b_itor) {
	  faceBoundary.push_back(*i_itor);
	}
      }
      if (debug > 1) {std::cout << "    boundary point " << *b_itor << std::endl;}
      this->buildFaces(dim-1, curSimplex, &faceBoundary, face);
      faces->push_back(face);
    }
  } else {
    if (debug > 1) {std::cout << "  Just set faces to boundary in 1d" << std::endl;}
    faces = boundary;
  }
  if (debug > 1) {
    for(PointArray::iterator f_itor = faces->begin(); f_itor != faces->end(); ++f_itor) {
      std::cout << "  face point " << *f_itor << std::endl;
    }
  }
  // We always create the toplevel, so we could short circuit somehow
  // Should not have to loop here since the meet of just 2 boundary elements is an element
  PointArray::iterator f_itor = faces->begin();
  point_type           start = *f_itor;
  f_itor++;
  point_type           next = *f_itor;
  Obj<sieve_type::supportSet> preElement = this->topology->nJoin(start, next, 1);
  
  if (preElement->size() > 0) {
    simplex = *preElement->begin();
    if (debug > 1) {std::cout << "  Found old simplex " << simplex << std::endl;}
  } else {
    int color = 0;
    
    simplex = point_type(0, (*(*curSimplex)[dim])++);
    for(PointArray::iterator f_itor = faces->begin(); f_itor != faces->end(); ++f_itor) {
      this->topology->addArrow(*f_itor, simplex, color++);
    }
    if (debug > 1) {std::cout << "  Added simplex " << simplex << " dim " << dim << std::endl;}
  }
  if (dim > 1) {
    delete faces;
  }
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::buildTopology"
// Build a topology from a connectivity description
//   (0, 0)            ... (0, numSimplices-1):  dim-dimensional simplices
//   (0, numSimplices) ... (0, numVertices):     vertices
// The other simplices are numbered as they are requested
void
ALE::PetscMesh::buildTopology(int numSimplices, 
			      int simplices[], 
			      int numVertices, 
			      bool interpolate,
			      int corners) {
  ALE_LOG_EVENT_BEGIN;
  if (this->commRank() != 0) {
    ALE_LOG_EVENT_END;
    return;
  }
  // Create a map from dimension to the current element number for that dimension
  std::map<int,int*> curElement = std::map<int,int*>();
  int                curSimplex = 0;
  int                curVertex  = numSimplices;
  int                newElement = numSimplices+numVertices;
  PointArray         boundary   = PointArray();
  
  if (corners < 0) corners = this->dim+1;
  curElement[0]         = &curVertex;
  curElement[this->dim] = &curSimplex;
  for(int d = 1; d < this->dim; d++) {
    curElement[d] = &newElement;
  }
  for(int s = 0; s < numSimplices; s++) {
    point_type simplex(0, s);
    
    // Build the simplex
    if (interpolate) {
      boundary.clear();
      for(int b = 0; b < corners; b++) {
	point_type vertex(0, simplices[s*corners+b]+numSimplices);
	
	if (debug > 1) {std::cout << "Adding boundary node " << vertex << std::endl;}
	boundary.push_back(vertex);
      }
      if (debug) {std::cout << "simplex " << s << " boundary size " << boundary.size() << std::endl;}
      
      if (corners != this->dim+1) {
	this->buildHexFaces(this->dim, &curElement, &boundary, simplex);
      } else {
	this->buildFaces(this->dim, &curElement, &boundary, simplex);
      }
    } else {
      for(int b = 0; b < corners; b++) {
	point_type p(0, simplices[s*corners+b]+numSimplices);
	
	this->topology->addArrow(p, simplex, b);
      }
    }
  }
  ALE_LOG_EVENT_END;
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::createVertBnd"
void
ALE::PetscMesh::createVertexBundle(int numSimplices, 
				   int simplices[], 
				   int elementOffset, 
				   int corners) {
  ALE_LOG_STAGE_BEGIN;
  Obj<bundle_type> vertexBundle = this->getBundle(0);
  Obj<sieve_type::traits::heightSequence> elements = this->topology->heightStratum(0);
  std::string orderName("element");
  int vertexOffset;
  
  ALE_LOG_EVENT_BEGIN;
  if (elementOffset) {
    vertexOffset = 0;
  } else {
    vertexOffset = numSimplices;
  }
  if (corners < 0) corners = this->dim+1;
  for(sieve_type::traits::heightSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); ++e_iter) {
    // setFiberDimensionByDepth() does not work here since we only want it to apply to the patch cone
    //   What we really need is the depthStratum relative to the patch
    Obj<PointArray> patch = PointArray();
    
    for(int b = 0; b < corners; b++) {
      patch->push_back(point_type(0, simplices[((*e_iter).index - elementOffset)*corners+b]+vertexOffset));
    }
    vertexBundle->setPatch(orderName, patch, *e_iter);
    for(PointArray::iterator p_iter = patch->begin(); p_iter != patch->end(); ++p_iter) {
      vertexBundle->setFiberDimension(orderName, *e_iter, *p_iter, 1);
    }
  }
  if (elements->size() == 0) {
    vertexBundle->setPatch(orderName, elements, bundle_type::patch_type());
  }
  ALE_LOG_EVENT_END;
  vertexBundle->orderPatches(orderName);
  ALE_LOG_STAGE_END;
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::createSerCoords"
void
ALE::PetscMesh::createSerialCoordinates(int embedDim, int numSimplices, double coords[]) {
  ALE_LOG_EVENT_BEGIN;
  patch_type patch;
  
  this->coordinates->setTopology(this->topology);
  this->coordinates->setPatch(this->topology->leaves(), patch);
  this->coordinates->setFiberDimensionByDepth(patch, 0, embedDim);
  this->coordinates->orderPatches();
  Obj<sieve_type::traits::depthSequence> vertices = this->topology->depthStratum(0);
  for(sieve_type::traits::depthSequence::iterator v_itor = vertices->begin(); v_itor != vertices->end(); v_itor++) {
    this->coordinates->update(patch, *v_itor, &coords[((*v_itor).index - numSimplices)*embedDim]);
  }
  Obj<bundle_type> vertexBundle = this->getBundle(0);
  Obj<sieve_type::traits::heightSequence> elements = this->topology->heightStratum(0);
  std::string orderName("element");
  
  for(sieve_type::traits::heightSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); e_iter++) {
    // setFiberDimensionByDepth() does not work here since we only want it to apply to the patch cone
    //   What we really need is the depthStratum relative to the patch
    Obj<bundle_type::order_type::coneSequence> cone = vertexBundle->getPatch(orderName, *e_iter);
    
    this->coordinates->setPatch(orderName, cone, *e_iter);
    for(bundle_type::order_type::coneSequence::iterator c_iter = cone->begin(); c_iter != cone->end(); ++c_iter) {
      this->coordinates->setFiberDimension(orderName, *e_iter, *c_iter, embedDim);
    }
  }
  if (elements->size() == 0) {
    this->coordinates->setPatch(orderName, elements, field_type::patch_type());
  }
  this->coordinates->orderPatches(orderName);
  ALE_LOG_EVENT_END;
}


#undef __FUNCT__
#define __FUNCT__ "Mesh::parVertReOrd"
// This is not right, we should not have to copy everything to the new order first
void
ALE::PetscMesh::createParallelVertexReorder(Obj<bundle_type> serialVertexBundle) {
  ALE_LOG_EVENT_BEGIN;
  Obj<bundle_type> vertexBundle = this->getBundle(0);
  std::string orderName("element");
  
  if (!this->commRank()) {
    Obj<bundle_type::order_type::baseSequence> patches = serialVertexBundle->getPatches(orderName);
    
    for(bundle_type::order_type::baseSequence::iterator e_iter = patches->begin(); e_iter != patches->end(); ++e_iter) {
      Obj<bundle_type::order_type::coneSequence> patch = serialVertexBundle->getPatch(orderName, *e_iter);
      
      vertexBundle->setPatch(orderName, patch, *e_iter);
      for(bundle_type::order_type::coneSequence::iterator p_iter = patch->begin(); p_iter != patch->end(); ++p_iter) {
	vertexBundle->setFiberDimension(orderName, *e_iter, *p_iter, 1);
      }
    }
  } else {
    Obj<sieve_type::traits::heightSequence> elements = this->topology->heightStratum(0);
    Obj<bundle_type::order_type> reorder = vertexBundle->__getOrder(orderName);
    
    for(sieve_type::traits::heightSequence::iterator e_iter = elements->begin(); e_iter != elements->end(); e_iter++) {
      reorder->addBasePoint(*e_iter);
    }
  }
  vertexBundle->orderPatches(orderName);
  vertexBundle->partitionOrder(orderName);
  ALE_LOG_EVENT_END;
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::populate"
// Create a serial mesh
void
ALE::PetscMesh::populate(int numSimplices, 
			 int simplices[], 
			 int numVertices, 
			 double coords[], 
			 bool interpolate, 
			 int corners) {
  this->topology->setStratification(false);
  this->buildTopology(numSimplices, simplices, numVertices, interpolate, corners);
  this->topology->stratify();
  this->topology->setStratification(true);
  this->createVertexBundle(numSimplices, simplices, 0, corners);
  this->createSerialCoordinates(this->dim, numSimplices, coords);
}

#undef __FUNCT__
#define __FUNCT__ "Mesh::populateBd"
void
ALE::PetscMesh::populateBd(int numSimplices, 
			   int simplices[], 
			   int numVertices, 
			   double coords[], 
			   bool interpolate, 
			   int corners) {
  this->topology->setStratification(false);
  this->buildTopology(numSimplices, simplices, numVertices, interpolate, corners);
  this->topology->stratify();
  this->topology->setStratification(true);
  this->createVertexBundle(numSimplices, simplices, 0, corners);
  this->createSerialCoordinates(this->dim+1, numSimplices, coords);
}

