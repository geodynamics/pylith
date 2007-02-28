pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt1.bc --pl3dscan.asciiOutputFile=pt1-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt2.bc --pl3dscan.asciiOutputFile=pt2-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt3.bc --pl3dscan.asciiOutputFile=pt3-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt4.bc --pl3dscan.asciiOutputFile=pt4-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt5.bc --pl3dscan.asciiOutputFile=pt5-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt6.bc --pl3dscan.asciiOutputFile=pt6-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt7.bc --pl3dscan.asciiOutputFile=pt7-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt8.bc --pl3dscan.asciiOutputFile=pt8-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt9.bc --pl3dscan.asciiOutputFile=pt9-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt10.bc --pl3dscan.asciiOutputFile=pt10-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt11.bc --pl3dscan.asciiOutputFile=pt11-bbar-comp.ascii
/bin/rm pt1.plot
pylith3dapp.py pt-bbar.cfg --pl3dscan.bcInputFile=pt12.bc --pl3dscan.asciiOutputFile=pt12-bbar-comp.ascii
/bin/rm pt1.plot
