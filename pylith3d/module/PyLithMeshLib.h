#ifdef __cplusplus
#define __PYX_EXTERN_C extern "C"
#else
#define __PYX_EXTERN_C extern
#endif
__PYX_EXTERN_C DL_IMPORT(PyTypeObject) PyMesh_Type;

struct PyMeshObject {
  PyObject_HEAD
  struct __pyx_vtabstruct_13PyLithMeshLib_Mesh *__pyx_vtab;
  Mesh mesh;
  Mat A;
  Vec rhs;
  Vec sol;
  char (*meshInputFile);
  char (*meshBcFile);
  int interpolateMesh;
  char (*partitioner);
  PyObject *_meshInputFile;
  PyObject *_meshBcFile;
  PyObject *_partitioner;
};
PyMODINIT_FUNC initPyLithMeshLib(void);
