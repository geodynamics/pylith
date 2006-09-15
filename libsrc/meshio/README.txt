Abstract

  public

    read(ALE::PetscMesh*)
    write(const ALE::PetscMesh&)

  protected, pure virtual

    _open(filename)
    _readMeshInfo(dimension)
    _writeMeshInfo(dimension)
    _readVertices(coordinates, numVertices, numDims)
    _writeVertices(coordinates, numVertices, numDims)
    _readElements(elements, numElements, numCorners)
    _writeElements(elements, numElements, numCorners)
    _readGroup()
    _writeGroup()
    _close()

  protected

    _addGroup(name, dimension, entries, numEntries)
    _getGroup(name, dimension, entries, numEntries)
