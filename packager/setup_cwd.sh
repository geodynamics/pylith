CWD=`pwd`
PATH=${CWD}/bin:/usr/bin:/bin
export PYTHONPATH=${CWD}/lib/python2.7/site-packages
export LD_LIBRARY_PATH=${CWD}/lib:${CWD}/lib64
unset PETSC_DIR PETSC_ARCH
