// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/CohesiveTopology.hh
 *
 * @brief C++ object to manage creation of cohesive cells.
 */

#if !defined(pylith_faults_cohesivetopology_hh)
#define pylith_faults_cohesivetopology_hh

#include "pylith/utils/sievetypes.hh" // USES PETSc Mesh

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class CohesiveTopology;
  } // faults
} // pylith

/// C++ object to manage creation of cohesive cells.
class pylith::faults::CohesiveTopology
{ // class CohesiveTopology
public :
  typedef std::set<Mesh::point_type>             PointSet;
  typedef std::vector<sieve_type::point_type>    PointArray;
  typedef std::pair<sieve_type::point_type, int> oPoint_type;
  typedef std::vector<oPoint_type>               oPointArray;
protected:
  template<typename Sieve, typename Renumbering>
  class ReplaceVisitor {
  public:
    typedef typename Sieve::point_type point_type;
  protected:
    Renumbering& renumbering;
    const int    size;
    int          i;
    const int    debug;
    point_type  *points;
    bool         mapped;
  public:
    ReplaceVisitor(Renumbering& r, const int size, const int debug = 0) : renumbering(r), size(size), i(0), debug(debug) {
      this->points = new point_type[this->size];
      this->mapped = false;
    };
    ~ReplaceVisitor() {delete [] this->points;};
    void visitPoint(const point_type& point) {
      if (i >= this->size) {throw ALE::Exception("Too many points for ReplaceVisitor");}
      if (this->renumbering.find(point) != this->renumbering.end()) {
        if (debug) std::cout << "    point " << this->renumbering[point] << std::endl;
        points[i] = this->renumbering[point];
        this->mapped = true;
      } else {
        if (debug) std::cout << "    point " << point << std::endl;
        points[i] = point;
      }
      ++i;
    };
    void visitArrow(const typename Sieve::arrow_type&) {};
  public:
    const point_type *getPoints() {return this->points;};
    bool mappedPoint() {return this->mapped;};
    void clear() {this->i = 0; this->mapped = false;};
  };
  template<typename Sieve>
  class ClassifyVisitor {
  public:
    typedef typename Sieve::point_type point_type;
  protected:
    const Sieve&     sieve;
    const PointSet&  replaceCells;
    const PointSet&  noReplaceCells;
    const point_type firstCohesiveCell;
    const int        faceSize;
    const int        debug;
    PointSet         vReplaceCells;
    PointSet         vNoReplaceCells;
    bool             modified;
    bool             setupMode;
    int              size;
    ALE::ISieveVisitor::PointRetriever<Sieve> pR;
  public:
    ClassifyVisitor(const Sieve& s, const PointSet& rC, const PointSet& nrC, const point_type& fC, const int fS, const int debug = 0) : sieve(s), replaceCells(rC), noReplaceCells(nrC), firstCohesiveCell(fC), faceSize(fS), debug(debug), modified(false), setupMode(true), size(0) {
      pR.setSize(s.getMaxConeSize());
    };
    ~ClassifyVisitor() {};
    void visitPoint(const point_type& point) {
      if (this->setupMode) {
        if (replaceCells.find(point)   != replaceCells.end())   vReplaceCells.insert(point);
        if (noReplaceCells.find(point) != noReplaceCells.end()) vNoReplaceCells.insert(point);
        if (point >= firstCohesiveCell) return;
        this->modified = true;
        this->size++;
        return;
      }
      bool classified = false;

      if (debug) {std::cout << "Checking neighbor " << point << std::endl;}
      if (vReplaceCells.find(point)   != vReplaceCells.end()) {
        if (debug) {std::cout << "  already in replaceCells" << std::endl;}
        return;
      }
      if (vNoReplaceCells.find(point) != vNoReplaceCells.end()) {
        if (debug) {std::cout << "  already in noReplaceCells" << std::endl;}
        return;
      }
      if (point >= firstCohesiveCell) {
        if (debug) {std::cout << "  already a cohesive cell" << std::endl;}
        return;
      }
      // If neighbor shares a face with anyone in replaceCells, then add
      for(PointSet::const_iterator c_iter = vReplaceCells.begin(); c_iter != vReplaceCells.end(); ++c_iter) {
        sieve.meet(*c_iter, point, pR);

        if (pR.getSize() == faceSize) {
          if (debug) {std::cout << "    Scheduling " << point << " for replacement" << std::endl;}
          vReplaceCells.insert(point);
          modified   = true;
          classified = true;
          pR.clear();
          break;
        }
        pR.clear();
      }
      if (classified) return;
      // It is unclear whether taking out the noReplace cells will speed this up
      for(PointSet::const_iterator c_iter = vNoReplaceCells.begin(); c_iter != vNoReplaceCells.end(); ++c_iter) {
        sieve.meet(*c_iter, point, pR);

        if (pR.getSize() == faceSize) {
          if (debug) {std::cout << "    Scheduling " << point << " for no replacement" << std::endl;}
          vNoReplaceCells.insert(point);
          modified   = true;
          classified = true;
          pR.clear();
          break;
        }
        pR.clear();
      }
    };
    void visitArrow(const typename Sieve::arrow_type&) {};
  public:
    const PointSet& getReplaceCells() const {return this->vReplaceCells;};
    const PointSet& getNoReplaceCells() const {return this->vNoReplaceCells;};
    const bool      getModified() const {return this->modified;};
    const int       getSize() const {return this->size;};
    void            setMode(const bool isSetup) {this->setupMode = isSetup;};
    void            reset() {this->modified = false;};
  };

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :
  /** Create the fault mesh.
   *
   * @param fault Finite-element mesh of fault (output)
   * @param mesh Finite-element mesh
   * @param faultVertices Vertices assocated with faces of cells defining 
   *   fault surface
   */
  static
  void createFault(Obj<SubMesh>& ifault,
                   Obj<ALE::Mesh>& faultBd,
                   const Obj<Mesh>& mesh,
                   const Obj<Mesh::int_section_type>& groupField);

  /** Create cohesive cells.
   *
   * @param fault Finite-element mesh of fault (output)
   * @param mesh Finite-element mesh
   * @param materialId Material id for cohesive elements.
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise
   */
  static
  void create(Obj<SubMesh>& ifault,
              const Obj<ALE::Mesh>& faultBd,
              const Obj<Mesh>& mesh,
              const Obj<Mesh::int_section_type>& groupField,
              const int materialId,
              const bool constraintCell = false,
              const bool flipFault = false);

  /** Create (distributed) fault mesh from cohesive cells.
   *
   * @param fault Finite-element mesh of fault (output).
   * @param cohesiveToFault Mapping of cohesive cell to fault mesh cell.
   * @param mesh Finite-element mesh.
   * @param materialId Material id for cohesive elements.
   * @param constraintCell True if creating cells constrained with 
   *   Lagrange multipliers that require extra vertices, false otherwise.
   */
  static
  void createParallel(ALE::Obj<SubMesh>* ifault,
		      std::map<Mesh::point_type, Mesh::point_type>* cohesiveToFault,
		      const ALE::Obj<Mesh>& mesh,
		      const int materialId,
		      const bool constraintCell =false);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :
  /** Get number of vertices on face.
   *
   * @param cell Finite-element cell
   * @param mesh Finite-element mesh
   *
   * @returns Number of vertices on cell face
   */
  static
  unsigned int _numFaceVertices(const Mesh::point_type& cell,
                                const ALE::Obj<Mesh>& mesh,
                                const int depth = -1);

  /** Determine a face orientation
   *    We should really have an interpolated mesh, instead of
   *    calculating this on the fly.
   *
   * @param cell Finite-element cell
   * @param mesh Finite-element mesh
   *
   * @returns True for positive orientation, otherwise false
   */
  static
  bool _faceOrientation(const Mesh::point_type& cell,
                        const ALE::Obj<Mesh>& mesh,
                        const int numCorners,
                        const int indices[],
                        const int oppositeVertex,
                        PointArray *origVertices,
                        PointArray *faceVertices);

  template<typename FaceType>
  static
  bool _getOrientedFace(const ALE::Obj<Mesh>& mesh,
                        const Mesh::point_type& cell,
                        FaceType face,
                        const int numCorners,
                        int indices[],
                        PointArray *origVertices,
                        PointArray *faceVertices);

  template<class InputPoints>
  static
  bool _compatibleOrientation(const ALE::Obj<Mesh>& mesh,
                              const Mesh::point_type& p,
                              const Mesh::point_type& q,
                              const int numFaultCorners,
                              const int faultFaceSize,
                              const int faultDepth,
                              const Obj<InputPoints>& points,
                              int indices[],
                              PointArray *origVertices,
                              PointArray *faceVertices,
                              PointArray *neighborVertices);

  static
  void _computeCensoredDepth(const ALE::Obj<Mesh::label_type>& depth,
                             const ALE::Obj<Mesh::sieve_type>& sieve,
                             const Mesh::point_type& firstCohesiveCell);

  static void classifyCells(const ALE::Obj<Mesh::sieve_type>& sieve,
                            const Mesh::point_type& vertex,
                            const int depth,
                            const int faceSize,
                            const Mesh::point_type& firstCohesiveCell,
                            PointSet& replaceCells,
                            PointSet& noReplaceCells,
                            const int debug);

  static void createFaultSieveFromVertices(const int dim,
                                           const int firstCell,
                                           const PointSet& faultVertices,
                                           const Obj<Mesh>& mesh,
                                           const Obj<ALE::Mesh::arrow_section_type>& orientation,
                                           const Obj<ALE::Mesh::sieve_type>& faultSieve);

  static void createFaultSieveFromFaces(const int dim,
                                        const int firstCell,
                                        const int numFaces,
                                        const int faultVertices[],
                                        const int faultCells[],
                                        const Obj<Mesh>& mesh,
                                        const Obj<ALE::Mesh::arrow_section_type>& orientation,
                                        const Obj<ALE::Mesh::sieve_type>& faultSieve);

  static void orientFaultSieve(const int dim,
                               const Obj<Mesh>& mesh,
                               const Obj<ALE::Mesh::arrow_section_type>& orientation,
                               const Obj<ALE::Mesh>& fault);
}; // class CohesiveTopology

#endif // pylith_faults_cohesivetopology_hh


// End of file 
