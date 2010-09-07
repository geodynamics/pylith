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

/**
 * @file libsrc/topology/UniformSectionDS.hh
 *
 * @brief Sieve section with uniform size set at runtime in contrast
 * to ALE::IUniformSection, which has a uniform size set at compile
 * time.
 */

#if !defined(pylith_topology_uniformsection_hh)
#define pylith_topology_uniformsection_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "IField.hh"

// UniformSectionDS -----------------------------------------------------
/// Sieve section with uniform size set at runtime.
/// All fibers are the same dimension
/// Note we can use a IConstantSection for this Atlas
/// Each point may have a different vector
/// Thus we need storage for values, and hence must implement completion
template<typename Point_, 
	 typename Value_, 
	 typename Alloc_ =ALE::malloc_allocator<Value_> >
class ALE::IUniformSectionDS : public ALE::ParallelObject {

// PUBLIC TYPEDEFS ////////////////////////////////////////////////////// 
public:
  typedef Point_ point_type;
  typedef Value_ value_type;
  typedef Alloc_ alloc_type;
  typedef typename alloc_type::template rebind<point_type>::other point_alloc_type;
  typedef IConstantSection<point_type, int, point_alloc_type> atlas_type;
  typedef typename atlas_type::chart_type chart_type;
  typedef point_type index_type;
  typedef struct { value_type* v;} fiber_type;
  typedef value_type* values_type;
  typedef typename alloc_type::template rebind<atlas_type>::other atlas_alloc_type;
  typedef typename atlas_alloc_type::pointer atlas_ptr;

// PUBLIC METHODS ///////////////////////////////////////////////////////
public:

  IUniformSectionDS(MPI_Comm comm, 
		    const int fiberDim,
		    const int debug =0);

  IUniformSectionDS(MPI_Comm comm, 
		    const int diberDim,
		    const point_type& min,
		    const point_type& max, 
		    const int debug =0);
  
  IUniformSectionDS(const Obj<atlas_type>& atlas,
		    const int fiberDim);
  
  /// Destructor.
  virtual ~IUniformSectionDS(void);

  value_type* getRawArray(const int size);

  bool hasPoint(const point_type& point) const;

  void checkDimension(const int& dim);

  const chart_type& getChart(void) const;

  void setChart(const chart_type& chart);

  bool resizeChart(const chart_type& chart);

  const Obj<atlas_type>& getAtlas(void) const;

  void setAtlas(const Obj<atlas_type>& atlas);
  void addPoint(const point_type& point);

  template<typename Points>
  void addPoint(const Obj<Points>& points);

  void copy(const Obj<IUniformSectionDS>& section);

  const value_type* getDefault(void) const;

  void setDefault(const value_type v[]);

  void clear(void);

  int getFiberDimension(const point_type& p) const;

  void setFiberDimension(const point_type& p,
			 int dim);

  template<typename Sequence>
  void setFiberDimension(const Obj<Sequence>& points,
			 int dim);

  void setFiberDimension(const std::set<point_type>& points,
			 int dim);

  void addFiberDimension(const point_type& p,
			 int dim);

  int size(void) const;
  
  int sizeWithBC(void) const;

  void allocatePoint(void);

  bool reallocatePoint(const chart_type& chart,
		       values_type *oldData =NULL);

  template<typename Iterator,
	   typename Extractor>
  bool reallocatePoint(const Iterator& begin,
		       const Iterator& end,
		       const Extractor& extractor);

  void zero(void);

  // Return a pointer to the entire contiguous storage array
  const values_type& restrictSpace(void) const;

  // Return only the values associated to this point, not its closure
  const value_type *restrictPoint(const point_type& p) const;

  // Update only the values associated to this point, not its closure
  void updatePoint(const point_type& p,
		   const value_type v[]);

  // Update only the values associated to this point, not its closure
  void updateAddPoint(const point_type& p,
		      const value_type v[]);

  void updatePointAll(const point_type& p, 
		      const value_type v[]);

  void view(const std::string& name,
	    MPI_Comm comm =MPI_COMM_NULL) const;

  int getNumSpaces(void) const;

  const std::vector<Obj<atlas_type> >& getSpaces(void);
  void addSpace(void);
  int getFiberDimension(const point_type& p,
			const int space) const;
  void setFiberDimension(const point_type& p,
			 int dim,
			 const int space);
  
  template<typename Sequence>
  void setFiberDimension(const Obj<Sequence>& points,
			 int dim,
			 const int space);

  // Return the total number of free dofs
  int size(const int space) const;

  template<typename OtherSection>
  void copyFibration(const Obj<OtherSection>& section);

  Obj<IUniformSectionDS> getFibration(const int space) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:
  int _fiberDim;
  Obj<atlas_type> _atlas;
  std::vector<Obj<atlas_type> > _spaces;
  values_type _array;
  fiber_type _emptyValue;
  alloc_type _allocator;

}; // class IUniformSectionDS

#include "UniformSectionDS.cc" // template definitions

#endif // uniformsectionds_hh

// End of file
