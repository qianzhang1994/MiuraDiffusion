"""
Main Procedure to conduct the parametric analysis
including CreaseArrayInput.py, RunABAQSU.py and odbWrite.py
The load is subjected to the part middle edge
"""
from CreaseArrayInput import CreaseArray
from RunABAQUS import RunABAQUS
from DataProcess import DataPostProcess
import os
import shutil
import pandas as pd
from math import *
import numpy as np
from scipy import stats

"""
Working directory settings, and copy the post process file (.py)
to the calculation directory
"""
PathNow = os.getcwd()
print("The Current program code path", PathNow)
Path = 'C:\ABAQUSTemp\MiuraArray\Simulation0210'
# shutil.rmtree(Path)
# os.mkdir(Path)
shutil.copyfile("OdbWriteStep.py", Path+"\OdbWriteStep.py")
os.chdir(Path)
PathNew = os.getcwd()
print("The Calculation analysis save path", PathNew)

"""
Input Modulus for parametric analysis
parameter Stiffness:    the crease stiffness per length
Parameter SectorAngle:  The plate geometry, [width, height, Sector Angle]
                        the crease length is the height, and  Sector angle between line 1/4 to 1/2, default = 90
Parameter InitialAngle: Define the initial rest angle of the single crease
Parameter dx:           mesh size
ArrayNumberX:           The count of MIURA unit in the X direction.
ArrayNumberY:           The count of MIURA unit in the Y direction.
Parameter foldingstate: 1 is unfolding and 0 is the folding process,default = 1

[No use Parameters in this Program]
Parameter Loading region:   for the ratio of load region to the edge set. 
                            It is the middle region. 0.5 means the middle 50% are loaded
CreaseRegion:               defined to expression the cut region. 
                            [0.2 0.5] expresses that there is no connection during the region [0.2 0.5]
"""

def Crease(jobnumber, Stiffness, height, SectorAngle, InitialAngle, width, thickness, ArrayNumberX, ArrayNumberY, CreaseNumber):
    jobname = "CreaseArray"+str(jobnumber)
    with open('readmejob.txt', "w") as f:
        f.write(jobname)
    dx = 5.0
    Model = CreaseArray(jobname, Stiffness, dx=dx, width=width, height=height,
                        SectorAngle=SectorAngle, InitialAngle=InitialAngle, thickness=thickness,
                        foldingState=1, Loadregion=1.0, Creaseregion=[0.8, 0.2], ArrayNumberX=ArrayNumberX,
                        ArrayNumberY=ArrayNumberY, CreaseNumber=CreaseNumber)
    Model.Inpwrite()
    # ABAQUSModel = RunABAQUS(jobname)
    # ABAQUSModel.Simulate()
    # for Step in range(1, 3):
    #     stepname = "Step-"+str(Step)
    #     with open('readmeStep.txt', "w") as f:
    #         f.write(stepname)
    #     os.system('abaqus cae noGUI=OdbWriteStep.py')
    #     DataModel = DataPostProcess(jobname, jobnumber, Stiffness, dx, InitialAngle)
    #     DataModel.Process()

"""
Main procedure for Crease Array analysis
ArrayNumberX is fixed to be 2
CreaseNumber = [1, MaxCreaseNumber], for the symmetry condition of ArrayY, CreaseNUmber = ArrayNumberY
foldingState: 1 is unfolding and 0 is the folding process,default = 1
if foldingstate = 0 and ArrayNumberY is odd (3, 5), U4 is positive (InitialAngle)
if foldingstate = 0 and ArrayNumberY is even (2, 4), U4 is negative (-InitialAngle)
if foldingstate = 1 and ArrayNumberY is odd (3, 5), U4 is negative (pi - InitialAngle)
if foldingstate = 1 and ArrayNumberY is even (2, 4), U4 is positive (pi - InitialAngle)
Local over-unfolding is not considered
"""

jobnumber = 209
# Stiff = [5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 5e-2, 5e-1]
# for Sector in [60, 85]:
#     for CreaseNumber in [1, 2, 3]:
for Sector in [85]:
    for CreaseNumber in range(9, 17):
        for Stiffness in [5e-4, ]:
            print("Jobnumber = ", jobnumber)

            # Crease Parameter
            print(Stiffness)

            # Geometry Parameters of the plates
            height = 100
            width = 100
            SectorAngle = Sector / 180 * pi
            thickness = 0.1

            # Array Parameters
            InitialAngle = 0.1 / 180 * pi
            ArrayNumberX = 2
            ArrayNumberY = CreaseNumber
            # ArrayNumberY = CreaseNumber
            NaxCreaseNumber = (ArrayNumberX - 1) * ArrayNumberY + ArrayNumberX * (ArrayNumberY - 1)
            print("NaxCreaseNumber = ", NaxCreaseNumber)

            # Loading Parameters and Boundary Conditions Modulus
            # CreaseNumber = 1
            foldingState = 1

            # Call the Main Calculation Modulus
            Crease(jobnumber, Stiffness, height, SectorAngle, InitialAngle, width, thickness,
                   ArrayNumberX, ArrayNumberY, CreaseNumber)
            jobnumber = jobnumber + 1



