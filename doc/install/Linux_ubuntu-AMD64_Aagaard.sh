#!/bin/bash                                                                                                              

#
if [ "$TOOLS_FORMAT" == "" ]; then
  export TOOLS_FORMAT=gcc-4.4.1_64
fi                                

# Store pre-tools info so we can reread this file after changing settings
if (($READTOOLS)); then                                                  
  PATH=${path_notools}                                                   
  export LD_LIBRARY_PATH=${ld_library_path_notools}                      
  export MANPATH=${manpath_notools}                                      
  export PYTHONPATH=${pythonpath_notools}                                
else                                                                     
  export path_notools=${PATH}                                            
  export ld_library_path_notools=${LD_LIBRARY_PATH}                      
  export manpath_notools=${MANPATH}                                      
  export pythonpath=${PYTHONPATH}                                        
  export READTOOLS=1                                                     
fi                                                                       

PYTHON_VERSION=2.6

export TOOLS_DIR=$HOME/tools/${TOOLS_FORMAT}
if [ "$LD_LIBRARY_PATH" == "" ]; then                
  export LD_LIBRARY_PATH=${TOOLS_DIR}/lib            
else                                                 
  export LD_LIBRARY_PATH=${TOOLS_DIR}/lib:${LD_LIBRARY_PATH}
fi                                                          
if [ "${PYTHON_PATH}" == "" ]; then                         
  export PYTHONPATH=${TOOLS_DIR}/lib/python${PYTHON_VERSION}/site-packages
else                                                                      
  export PYTHONPATH=${TOOLS_DIR}/lib/python${PYTHON_VERSION}/site-packages:${PYTHONPATH}
fi
PATH=${TOOLS_DIR}/bin:${PATH}

# PETSC
export PETSC_DIR=$HOME/src/petsc-dev
export PETSC_ARCH=${TOOLS_FORMAT}
