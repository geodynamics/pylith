pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test1.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test1-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test2.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test2-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test3.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test3-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test4.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test4-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test5.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test5-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test6.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test6-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test7.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test7-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test8.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test8-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test9.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test9-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test10.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test10-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test11.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test11-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
pycrust3dapp.py --pc3dscan.fileRoot=pt --pc3dscan.keywordEqualsValueFile=pt-2pt.keyval --pc3dscan.bcInputFile=pt-test12.bc --pc3dscan.timeStepInputFile=pt-nobbar.time --pc3dscan.asciiOutputFile=pt-test12-2pt-nobbar-incomp.ascii --pc3dscan.materialPropertiesInputFile=pt-incomp.prop
/bin/rm pt.plot
