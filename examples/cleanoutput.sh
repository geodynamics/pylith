#!/bin/bash

for suffix in h5 xmf vtk; do
  rm `find . -name "*.$suffix"`
done

exit 0
