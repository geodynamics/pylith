pylith3dapp.py --bcInputFile=pt1.bc --asciiOutputFile=pt1-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt2.bc --asciiOutputFile=pt2-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt3.bc --asciiOutputFile=pt3-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt4.bc --asciiOutputFile=pt4-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt5.bc --asciiOutputFile=pt5-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt6.bc --asciiOutputFile=pt6-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt7.bc --asciiOutputFile=pt7-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt8.bc --asciiOutputFile=pt8-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt9.bc --asciiOutputFile=pt9-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt10.bc --asciiOutputFile=pt10-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt11.bc --asciiOutputFile=pt11-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py --bcInputFile=pt12.bc --asciiOutputFile=pt12-full-comp.ascii
/bin/rm pt1.plot
