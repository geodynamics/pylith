pylith3dapp.py pt-bbar.cfg --bcInputFile=pt1.bc --asciiOutputFile=pt1-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt2.bc --asciiOutputFile=pt2-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt3.bc --asciiOutputFile=pt3-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt4.bc --asciiOutputFile=pt4-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt5.bc --asciiOutputFile=pt5-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt6.bc --asciiOutputFile=pt6-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt7.bc --asciiOutputFile=pt7-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt8.bc --asciiOutputFile=pt8-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt9.bc --asciiOutputFile=pt9-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt10.bc --asciiOutputFile=pt10-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt11.bc --asciiOutputFile=pt11-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --bcInputFile=pt12.bc --asciiOutputFile=pt12-bbar-comp.ascii
/bin/rm pt1.plot
