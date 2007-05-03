pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt1.bc --asciiOutputFile=pt1-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt2.bc --asciiOutputFile=pt2-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt3.bc --asciiOutputFile=pt3-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt4.bc --asciiOutputFile=pt4-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt5.bc --asciiOutputFile=pt5-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt6.bc --asciiOutputFile=pt6-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt7.bc --asciiOutputFile=pt7-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt8.bc --asciiOutputFile=pt8-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt9.bc --asciiOutputFile=pt9-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt10.bc --asciiOutputFile=pt10-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt11.bc --asciiOutputFile=pt11-red-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --fileRoot=pt1 --keywordEqualsValueFile=pt-red.keyval --bcInputFile=pt12.bc --asciiOutputFile=pt12-red-comp.ascii
/bin/rm pt1.plot
