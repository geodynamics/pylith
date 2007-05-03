pylith3dapp.py pt-full.cfg --bcInputFile=pt1.bc --asciiOutputFile=pt1-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt2.bc --asciiOutputFile=pt2-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt3.bc --asciiOutputFile=pt3-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt4.bc --asciiOutputFile=pt4-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt5.bc --asciiOutputFile=pt5-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt6.bc --asciiOutputFile=pt6-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt7.bc --asciiOutputFile=pt7-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt8.bc --asciiOutputFile=pt8-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt9.bc --asciiOutputFile=pt9-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt10.bc --asciiOutputFile=pt10-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt11.bc --asciiOutputFile=pt11-full-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-full.cfg --bcInputFile=pt12.bc --asciiOutputFile=pt12-full-comp.ascii
/bin/rm pt1.plot
