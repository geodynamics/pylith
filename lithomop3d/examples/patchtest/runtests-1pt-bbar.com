pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test1.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test1-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test2.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test2-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test3.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test3-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test4.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test4-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test5.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test5-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test6.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test6-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test7.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test7-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test8.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test8-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test9.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test9-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test10.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test10-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test11.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test11-1pt-bbar.ascii
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-1pt.keyval --pc3dscan.bcInputFile=pt-test12.bc --pc3dscan.timeStepInputFile=pt-bbar.time --pc3dscan.asciiOutputFile=pt-test12-1pt-bbar.ascii
/bin/rm pt.plot
