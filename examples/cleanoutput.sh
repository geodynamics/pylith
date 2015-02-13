#!/bin/bash

for suffix in h5 xmf vtk; do
  rm `find . -name "*.$suffix"`
done

rm `find . -name progress.txt`

exit 0
