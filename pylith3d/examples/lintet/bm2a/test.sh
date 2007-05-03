#!/bin/sh

pylith3dapp.py -fileRoot=lintet/bm1a/bm1a -bcInputFile=lintet/bm2a/bm2a.bc -keywordEqualsValueFile=lintet/bm2a/bm2a.keyval -timeStepInputFile=lintet/bm2a/bm2a.time -asciiOutputFile=lintet/bm2a/bm2a.ascii -ucdOutputRoot=lintet/bm2a/bm2a

# end of file
