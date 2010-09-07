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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::IUniformSectionDS(MPI_Comm comm,
					  const int fiberDim,
					  const int debug) : 
  ParallelObject(comm, debug)
{ // constructor
  if (fiberDim <= 0) {
    std::ostringstream msg;
    msg << "Fiber dimension '" << fiberDim << "' for section must be positive.";
    throw ALE::Exception(msg.str().c_str());
  } // if
  _fiberDim = fiberDim;

  atlas_ptr pAtlas = atlas_alloc_type(this->_allocator).allocate(1);
  atlas_alloc_type(this->_allocator).construct(pAtlas, atlas_type(comm, debug));
  this->_atlas = Obj<atlas_type>(pAtlas, sizeof(atlas_type));
  this->_array = NULL;
  this->_emptyValue.v = new value_type[fiberDim];
  for(int i = 0; i < fiberDim; ++i)
    this->_emptyValue.v[i] = value_type();
} // constructor

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::IUniformSectionDS(MPI_Comm comm,
					  const int fiberDim,
					  const point_type& min,
					  const point_type& max,
					  const int debug) :
  ParallelObject(comm, debug)
{ // constructor
  if (fiberDim <= 0) {
    std::ostringstream msg;
    msg << "Fiber dimension '" << fiberDim << "' for section must be positive.";
    throw ALE::Exception(msg.str().c_str());
  } // if
  _fiberDim = fiberDim;

  atlas_ptr pAtlas = atlas_alloc_type(this->_allocator).allocate(1);
  atlas_alloc_type(this->_allocator).construct(pAtlas, 
					       atlas_type(comm, min, max, fiberDim, debug));
  this->_atlas = Obj<atlas_type>(pAtlas, sizeof(atlas_type));
  this->_array = NULL;
  this->_emptyValue.v = new value_type[fiberDim];
  for(int i = 0; i < fiberDim; ++i)
    this->_emptyValue.v[i] = value_type();
} // constructor

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::IUniformSectionDS(const Obj<atlas_type>& atlas,
					const int fiberDim) :
  ParallelObject(atlas->comm(), atlas->debug()), _atlas(atlas)
{ // constructor
  if (fiberDim <= 0) {
    std::ostringstream msg;
    msg << "Fiber dimension '" << fiberDim << "' for section must be positive.";
    throw ALE::Exception(msg.str().c_str());
  } // if
  _fiberDim = fiberDim;

  this->_atlas->update(*this->_atlas->getChart().begin(), &fiberDim);
  this->_array = NULL;
  this->_emptyValue.v = new value_type[fiberDim];
  for(int i = 0; i < fiberDim; ++i)
    this->_emptyValue.v[i] = value_type();
} // constructor

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::~IUniformSectionDS(void)
{ // destructor
  delete[] this->_emptyValue.v; this->_emptyValue.v = NULL;

  if (this->_array) {
    const int chartEnd = this->getChart().max()*_fiberDim;
    for(int i = this->getChart().min()*_fiberDim;
	i < chartEnd;
	++i)
      this->_allocator.destroy(this->_array+i);
    this->_array += this->getChart().min()*_fiberDim;
    this->_allocator.deallocate(this->_array, this->sizeWithBC());
    this->_array = NULL;
    this->_atlas = NULL;
  } // if
} // destructor

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::value_type*
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getRawArray(const int size)
{ // getRawArray
  static value_type* array = NULL;
  static int maxSize = 0;

  if (size > maxSize) {
    const value_type dummy(0);

    if (array) {
      for(int i = 0; i < maxSize; ++i)
	this->_allocator.destroy(array+i);
      this->_allocator.deallocate(array, maxSize);
    } // if
    maxSize = size;
    array = this->_allocator.allocate(maxSize);
    for (int i = 0; i < maxSize; ++i)
      this->_allocator.construct(array+i, dummy);
  }
  return array;
} // getRawArray

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
bool
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::hasPoint(const point_type& point) const
{ // hasPoint
  return this->_atlas->hasPoint(point);
} // hasPoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::checkDimension(const int& dim)
{ // checkDimension
  if (dim != _fiberDim) {
    ostringstream msg;
    msg << "Invalid fiber dimension '" << dim << "' must be '" << _fiberDim
	<< "'." << std::endl;
    throw ALE::Exception(msg.str().c_str());
  } // if
} // checkDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::chart_type&
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getChart(void) const
{ // getChart
  return this->_atlas->getChart();
} // getChart

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setChart(const chart_type& chart)
{ // setChart
  this->_atlas->setChart(chart);
  int dim = _fiberDim;
  this->_atlas->updatePoint(*this->getChart().begin(), &dim);
} // setChart

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
bool
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::resizeChart(const chart_type& chart)
{ // resizeChart
  if ((chart.min() >= this->getChart().min()) &&
      (chart.max() <= this->getChart().max()))
    return false;
  this->setChart(chart);
  return true;
} // resizeChart

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const ALE::Obj<typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::atlas_type>&
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getAtlas(void) const
{ // getAtlas
  return this->_atlas;
} // getAtlas

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setAtlas(const Obj<atlas_type>& atlas)
{ // setAtlas
  this->_atlas = atlas;
} // setAtlas

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::addPoint(const point_type& point)
{ // addPoint
  this->setFiberDimension(point, _fiberDim);
} // addPoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
template<typename Points>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::addPoint(const Obj<Points>& points) {
  const typename Points::const_iterator pointsEnd = points->end();
  for(typename Points::iterator p_iter=points->begin();
      p_iter != pointsEnd;
      ++p_iter)
    this->setFiberDimension(*p_iter, _fiberDim);
} // addPoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::copy(const Obj<IUniformSectionDS>& section)
{ // copy
  this->getAtlas()->copy(section->getAtlas());
  const chart_type& chart = section->getChart();
  
  const typename chart_type::const_iterator chartEnd = chart.end();
  for(typename chart_type::const_iterator c_iter=chart.begin();
      c_iter != chartEnd;
      ++c_iter)
    this->updatePoint(*c_iter, section->restrictPoint(*c_iter));
} // copy

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::value_type*
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getDefault(void) const
{ // getDefault
  return this->_emptyValue.v;
} // getDefault

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setDefault(const value_type v[])
{ // setDefault
  for (int i = 0; i < _fiberDim; ++i)
    this->_emptyValue.v[i] = v[i];
} // setDefault

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::clear(void)
{ // clear
  this->zero();
  this->_atlas->clear();
} // clear

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getFiberDimension(const point_type& p) const
{ // getFiberDimension
  return this->_atlas->restrictPoint(p)[0];
} // getFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void 
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setFiberDimension(const point_type& p,
					  int dim)
{ // setFiberDimension
  this->checkDimension(dim);
  this->_atlas->addPoint(p);
  this->_atlas->updatePoint(p, &dim);
} // setFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
template<typename Sequence>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setFiberDimension(const Obj<Sequence>& points,
					  int dim)
{ // setFiberDimension
  const typename Sequence::const_iterator pointsEnd = points->end();
  for (typename Sequence::iterator p_iter=points->begin();
       p_iter != pointsEnd;
       ++p_iter)
    this->setFiberDimension(*p_iter, dim);
} // setFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setFiberDimension(const std::set<point_type>& points,
					  int dim)
{ // setFiberDimension
  const typename std::set<point_type>::const_iterator pointsEnd = points.end();
  for (typename std::set<point_type>::iterator p_iter=points.begin();
       p_iter != pointsEnd;
       ++p_iter)
    this->setFiberDimension(*p_iter, dim);
} // setFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::addFiberDimension(const point_type& p,
					  int dim)
{ // addFiberDimension
  if (this->hasPoint(p)) {
    ostringstream msg;
    msg << "Invalid addition to fiber dimension " << dim
	<< " cannot exceed " << _fiberDim << std::endl;
    throw ALE::Exception(msg.str().c_str());
  } else
    this->setFiberDimension(p, dim);
} // addFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::size(void) const
{ // size
  return this->_atlas->getChart().size()*_fiberDim;
} // size

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::sizeWithBC(void) const
{ // sizeWithBC
  return this->size();
} // sizeWithBC

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::allocatePoint(void)
{ // allocatePoint
  this->_array = this->_allocator.allocate(this->sizeWithBC());
  this->_array -= this->getChart().min()*_fiberDim;
  const index_type chartEnd = this->getChart().max()*_fiberDim;
  for(index_type i=this->getChart().min()*_fiberDim;
      i < chartEnd;
      ++i)
    this->_allocator.construct(this->_array+i, this->_emptyValue.v[0]);
} // allocatePoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
bool
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::reallocatePoint(const chart_type& chart,
					values_type* oldData)
{ // reallocatePoint
  const chart_type  oldChart = this->getChart();
  const int         oldSize  = this->sizeWithBC();
  values_type       oldArray = this->_array;
  if (!this->resizeChart(chart))
    return false;
  const int         size     = this->sizeWithBC();
  
  this->_array = this->_allocator.allocate(size);
  this->_array -= this->getChart().min()*_fiberDim;
  const index_type chartEnd = this->getChart().max()*_fiberDim;
  for(index_type i=this->getChart().min()*_fiberDim;
      i < chartEnd;
      ++i)
    this->_allocator.construct(this->_array+i, this->_emptyValue.v[0]);
  
  const index_type oldChartEnd = oldChart.max()*_fiberDim;
  for(index_type i=oldChart.min()*_fiberDim; i < oldChartEnd; ++i)
    this->_array[i] = oldArray[i];
  if (!oldData) {
    for(index_type i=oldChart.min()*_fiberDim;
	i < oldChartEnd;
	++i)
      this->_allocator.destroy(oldArray+i);
    oldArray += this->getChart().min()*_fiberDim;
    this->_allocator.deallocate(oldArray, oldSize);
    ///std::cout << "Freed IUniformSection data" << std::endl;
  } else {
    ///std::cout << "Did not free IUniformSection data" << std::endl;
    *oldData = oldArray;
  } // if/else
  return true;
} // reallocatePoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
template<typename Iterator, typename Extractor>
bool
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::reallocatePoint(const Iterator& begin,
					const Iterator& end,
					const Extractor& extractor)
{ // reallocatePoint
  point_type min = this->getChart().min();
  point_type max = this->getChart().max()-1;

  for(Iterator p_iter = begin; p_iter != end; ++p_iter) {
    min = std::min(extractor(*p_iter), min);
    max = std::max(extractor(*p_iter), max);
  } // for
  return reallocatePoint(chart_type(min, max+1));
} // reallocatePoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::zero(void)
{ // zero
  memset(this->_array+(this->getChart().min()*_fiberDim), 0,
	 this->sizeWithBC()*sizeof(value_type));
} // zero

// ----------------------------------------------------------------------
// Return a pointer to the entire contiguous storage array
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::values_type& 
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::restrictSpace(void) const
{ // restrictSpace
  return this->_array;
} // restrictSpace

// ----------------------------------------------------------------------
// Return only the values associated to this point, not its closure
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::value_type*
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::restrictPoint(const point_type& p) const
{ // restrictPoint
  if (!this->hasPoint(p))
    return this->_emptyValue.v;
  return &this->_array[p*_fiberDim];
} // restrictPoint

// ----------------------------------------------------------------------
// Update only the values associated to this point, not its closure
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::updatePoint(const point_type& p,
				    const value_type v[])
{ // updatePoint
  for(int i = 0, idx = p*_fiberDim; i < _fiberDim; ++i, ++idx)
    this->_array[idx] = v[i];
} // updatePoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
// Update only the values associated to this point, not its closure
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::updateAddPoint(const point_type& p,
				       const value_type v[])
{ // updateAddPoint
  for(int i = 0, idx = p*_fiberDim; i < _fiberDim; ++i, ++idx)
    this->_array[idx] += v[i];
} // updateAddPoint

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void 
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::updatePointAll(const point_type& p,
				       const value_type v[])
{ // updatePointAll
  this->updatePoint(p, v);
} // updatePointAll

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::view(const std::string& name,
			     MPI_Comm comm) const
{ // view
  ostringstream txt;
  int rank;
  
  if (comm == MPI_COMM_NULL) {
    comm = this->comm();
    rank = this->commRank();
  } else {
    MPI_Comm_rank(comm, &rank);
  }
  if (name == "") {
    if(rank == 0) {
      txt << "viewing an IUniformSection" << std::endl;
    }
  } else {
    if(rank == 0) {
      txt << "viewing IUniformSection '" << name << "'" << std::endl;
    }
  }
  const typename atlas_type::chart_type& chart = this->_atlas->getChart();
  values_type                            array = this->_array;
  
  for(typename atlas_type::chart_type::const_iterator p_iter = chart.begin(); p_iter != chart.end(); ++p_iter) {
    const int idx = (*p_iter)*_fiberDim;
    
    if (_fiberDim != 0) {
      txt << "[" << this->commRank() << "]:   " << *p_iter << " dim " << _fiberDim << "  ";
      for(int i = 0; i < _fiberDim; i++) {
	txt << " " << array[idx+i];
      }
      txt << std::endl;
    }
  }
  if (chart.size() == 0) {
    txt << "[" << this->commRank() << "]: empty" << std::endl;
  }
  PetscSynchronizedPrintf(comm, txt.str().c_str());
  PetscSynchronizedFlush(comm);
} // view

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getNumSpaces(void) const
{ // getNumSpaces
  return this->_spaces.size();
} // getNumSpaces

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
const std::vector<ALE::Obj<typename ALE::IUniformSectionDS<point_type, value_type, alloc_type>::atlas_type> >&
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getSpaces(void)
{ // getSpaces
  return this->_spaces;
} // getSpaces

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::addSpace(void)
{ // addSpace
  Obj<atlas_type> space = new atlas_type(this->comm(), this->debug());
  this->_spaces.push_back(space);
} // addSpace

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int 
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getFiberDimension(
						    const point_type& p,
						    const int space) const
{ // getFiberDimension
  return *this->_spaces[space]->restrictPoint(p);
} // getFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setFiberDimension(
					    const point_type& p,
					    int dim,
					    const int space)
{ // setFiberDimension
  const index_type idx = -1;
  this->_spaces[space]->addPoint(p);
  this->_spaces[space]->updatePoint(p, &idx);
} // setFiberDimension

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
template<typename Sequence>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::setFiberDimension(
					const Obj<Sequence>& points,
					int dim,
					const int space)
{ // setFiberDimension
  const typename Sequence::const_iterator pointsEnd = points->end();
  for(typename Sequence::iterator p_iter = points->begin();
      p_iter != pointsEnd;
      ++p_iter)
    this->setFiberDimension(*p_iter, dim, space);
} // setFiberDimension

// ----------------------------------------------------------------------
// Return the total number of free dofs
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
int
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::size(const int space) const
{ // size
  const chart_type& points = this->getChart();
  int size   = 0;
  
  const typename chart_type::const_iterator pointsEnd = points.end();
  for (typename chart_type::const_iterator p_iter = points.begin();
       p_iter != pointsEnd;
       ++p_iter)
    size += this->getConstrainedFiberDimension(*p_iter, space);

  return size;
} // size

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
template<typename OtherSection>
void
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::copyFibration(const Obj<OtherSection>& section)
{ // copyFibration
  const std::vector<Obj<atlas_type> >& spaces = section->getSpaces();
  
  this->_spaces.clear();
  const typename std::vector<Obj<atlas_type> >::const_iteraor spacesEnd = spaces.end();
  for(typename std::vector<Obj<atlas_type> >::const_iterator s_iter=spaces.begin();
      s_iter != spacesEnd;
      ++s_iter)
    this->_spaces.push_back(*s_iter);
} // copyFibration

// ----------------------------------------------------------------------
template<typename point_type, 
	 typename value_type, 
	 typename alloc_type>
ALE::Obj<ALE::IUniformSectionDS<point_type, value_type, alloc_type> >
ALE::IUniformSectionDS<point_type, value_type, alloc_type>::getFibration(const int space) const {
  const chart_type& chart = this->getChart();

  const int localFiberDim = this->getFiberDimension(*chart.begin(), space);
  int fiberDim = 0;
  MPI_Allreduce((void *) &localFiberDim, (void *) &fiberDim, 1, 
		MPI_INT, MPI_MAX, this->comm());
  assert(fiberDim > 0);
  Obj<IUniformSectionDS> field = new IUniformSectionDS(this->comm(),
						       fiberDim,
						       *chart.begin(),
						       *chart.end(),
						       this->debug());
  
  // Copy sizes
  const typename chart_type::const_iterator chartEnd = chart.end();
  for(typename chart_type::const_iterator c_iter=chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    const int fDim = this->getFiberDimension(*c_iter, space);
    if (fDim)
      field->setFiberDimension(*c_iter, fDim);
  } // for
  
  // Copy values
  const chart_type& newChart = field->getChart();
  const typename chart_type::const_iterator newChartEnd = newChart.end();
  for(typename chart_type::const_iterator c_iter = newChart.begin();
      c_iter != newChartEnd;
      ++c_iter)
    field->updatePoint(*c_iter, this->restrictPoint(*c_iter));

  return field;
} // getFibration


// End of file
