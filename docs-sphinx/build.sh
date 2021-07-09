#!/bin/bash

formats=("html")
if [ $# -gt 0 ]; then
    formats="$@"
fi

SRC_DIR=.
BUILD_DIR=_build

for format in ${formats[*]}; do
    sphinx-build -b ${format} ${SRC_DIR} ${BUILD_DIR}/${format}
done

exit 0
