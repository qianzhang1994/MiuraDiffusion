from abaqus import *
from odbAccess import *

"""
Read the job Information
and the Step Information
Set-1 corresponds to the unfolding to the plane
Set-2 corresponds to the external tension of the plane
"""
f = open('readmejob.txt', 'r')
req = f.readline()
f.close()
path = req + '.odb'
an_odb_object = openOdb(path=path)

f = open('readmeStep.txt', 'r')
step = f.readline()
f.close()
StepNumber = step

f = open('Output_Connection.txt', 'r')
connection = f.readline()
f.close()
ConnectionNumber = int(connection)


"""
The ALLSE of the plate and the total systems
as PlateASSE and TotalALLSE
"""
Platedata = an_odb_object.steps[StepNumber].historyRegions["ElementSet SET-PLATE"].historyOutputs['ALLSE'].data
Totaldata = an_odb_object.steps[StepNumber].historyRegions["Assembly ASSEMBLY"].historyOutputs['ALLSE'].data
Framenumber = len(Platedata)
DataFile = open('Data_PlateALLSE.txt', 'w')
DataFile.write('time,PlateALLSE\n')
for time, Data in Platedata:
    DataFile.write('%10.4E, %10.4E\n' % (time, Data))
DataFile.close()

DataFile = open('Data_TotalALLSEl.txt', 'w')
DataFile.write('TotalALLSE\n')
for time, Data in Totaldata:
    DataFile.write('%10.4E\n' % (Data))
DataFile.close()

"""
Record all the output result documents name in a file
Resultfile_Label
"""
ResultLabel = open('Resultfile_Label.txt','w')
ResultLabel.write("Data_PlateALLSE\n")
ResultLabel.write("Data_TotalALLSEl\n")
ResultLabel.close()



"""
OutputX for the Crease coordinates in field output
"""
for num in range(ConnectionNumber):
    print(num)
    countX = len(open("OutputX" + str(num + 1) + ".txt").readlines())
    OutPutXList = []
    for line in open("OutputX" + str(num + 1) + ".txt").readlines():
        OutPutXList.append(line.strip("\n"))
    for i in range(countX):
        if (i == 0) | (i == countX - 1):
            ResultLabel = open('Resultfile_Label.txt', 'a+')
            ResultLabel.write("Data_CoordX" + str(num * countX + i + 1) + "\n")
            ResultLabel.close()

            DataFile = open("Data_CoordX" + str(num * countX + i + 1) + ".txt", 'w')
            NodeLabel = OutPutXList[i]
            NodeLabel = int(NodeLabel)
            line = "XCoordX." + str(NodeLabel) + ",XCoordY." + str(NodeLabel) + ",XCoordZ." + str(NodeLabel)
            DataFile.write(line)
            DataFile.write('\n')
            for k in range(Framenumber):
                Creasedata = an_odb_object.steps[StepNumber].frames[k].fieldOutputs['COORD'].values[int(NodeLabel) - 1].data
                line = str(Creasedata[0]) + "," + str(Creasedata[1]) + "," + str(Creasedata[2])
                DataFile.write(line)
                DataFile.write('\n')
        DataFile.close()
#
# """
# OutputY: Connection Pair moment output (CM1) in history output
# """
#
# for num in range(ConnectionNumber):
#     countY = len(open("OutputLY" + str(num + 1) + ".txt").readlines())
#     OutPutYList = []
#     for line in open("OutputLY" + str(num + 1) + ".txt").readlines():
#         OutPutYList.append(line.strip("\n"))
#     for i in range(countY):
#         ResultLabel = open('Resultfile_Label.txt', 'a+')
#         ResultLabel.write("Data_URL1_" + str(num * countY + i + 1) + "\n")
#         ResultLabel.close()
#
#         DataFile = open("Data_URL1_" + str(num * countY + i + 1) + ".txt", 'w')
#         NodeLabel = OutPutYList[i]
#         NodeLabel = int(NodeLabel)
#         line = "URL1." + str(NodeLabel)
#         DataFile.write(line)
#         DataFile.write('\n')
#         for k in range(Framenumber):
#             dtm = an_odb_object.rootAssembly.datumCsyses['ASSEMBLY_CSYS-CON'+str(num+1)]
#             print(dtm)
#             fieldout = an_odb_object.steps[StepNumber].frames[k].fieldOutputs["UR"].getTransformedField(datumCsys=dtm)
#             Creasedata = fieldout.values[int(NodeLabel) - 1].data
#             line = str(Creasedata[0])
#             DataFile.write(line)
#             DataFile.write('\n')
#     DataFile.close()
#
#
# for num in range(ConnectionNumber):
#     countY = len(open("OutputRY" + str(num + 1) + ".txt").readlines())
#     OutPutYList = []
#     for line in open("OutputRY" + str(num + 1) + ".txt").readlines():
#         OutPutYList.append(line.strip("\n"))
#     for i in range(countY):
#         ResultLabel = open('Resultfile_Label.txt', 'a+')
#         ResultLabel.write("Data_URR1_" + str(num * countY + i + 1) + "\n")
#         ResultLabel.close()
#
#         DataFile = open("Data_URR1_" + str(num * countY + i + 1) + ".txt", 'w')
#         NodeLabel = OutPutYList[i]
#         NodeLabel = int(NodeLabel)
#         line = "URR1." + str(NodeLabel)
#         DataFile.write(line)
#         DataFile.write('\n')
#         for k in range(Framenumber):
#             dtm = an_odb_object.rootAssembly.datumCsyses['ASSEMBLY_CSYS-CON'+str(num+1)]
#             print(dtm)
#             fieldout = an_odb_object.steps[StepNumber].frames[k].fieldOutputs["UR"].getTransformedField(datumCsys=dtm)
#             Creasedata = fieldout.values[int(NodeLabel) - 1].data
#             line = str(Creasedata[0])
#             DataFile.write(line)
#             DataFile.write('\n')
#     DataFile.close()
#
#
# for num in range(ConnectionNumber):
#     countY = len(open("OutputRY" + str(num + 1) + ".txt").readlines())
#     OutPutYList = []
#     for line in open("OutputRY" + str(num + 1) + ".txt").readlines():
#         OutPutYList.append(line.strip("\n"))
#     for i in range(countY):
#         ResultLabel = open('Resultfile_Label.txt', 'a+')
#         ResultLabel.write("Data_CUR1_" + str(num * countY + i + 1) + "\n")
#         ResultLabel.close()
#
#         DataFile = open("Data_CUR1_" + str(num * countY + i + 1) + ".txt", 'w')
#         NodeLabel = OutPutYList[i]
#         NodeLabel = int(NodeLabel)
#         line = "CUR1." + str(NodeLabel)
#         DataFile.write(line)
#         DataFile.write('\n')
#         for k in range(Framenumber):
#             # dtm = an_odb_object.rootAssembly.datumCsyses['ASSEMBLY_CSYS-CON'+str(num+1)]
#             # print(dtm)
#             fieldout = an_odb_object.steps[StepNumber].frames[k].fieldOutputs["CUR"]
#             Creasedata = fieldout.values[int(NodeLabel) - 1].data
#             line = str(Creasedata[0])
#             DataFile.write(line)
#             DataFile.write('\n')
#     DataFile.close()
#

#
"""
OutputY: Connection Pair Relative Displacement output (CUR1) in history output
"""
for num in range(ConnectionNumber):
    countY = len(open("OutputCrease" + str(num + 1) + ".txt").readlines())
    for i in range(countY):
        ResultLabel = open('Resultfile_Label.txt', 'a+')
        ResultLabel.write("Data_CreasePair" + str(num * countY + i + 1) + "\n")
        ResultLabel.close()

        DataFile = open("Data_CreasePair" + str(num * countY + i + 1) + ".txt", 'w')
        CUR1data = an_odb_object.steps[StepNumber].historyRegions[
            "Element ASSEMBLY." + str(num * countY + i + 1)].historyOutputs['CUR1'].data
        DataFile.write("Element ASSEMBLY." + str(num * countY + i + 1)+"\n")
        for time, Data in CUR1data:
            DataFile.write('%10.4E\n' % Data)
        DataFile.close()

