"""
Data Processing Modulus
Combine all result in Pandas.DataFrame
df: columns
time/plateAllSE/TotalALLSE/XCOORD/YCOORD/Conection Pair Moment
XCOORD: the coordinates of the nodes in middle plate 0 along the width directin
YCOORD: the coordinate of the nodes in the connection edge of plate 0
XCOORD and YCOORD are calculated to determine the configuration.
counX and CountY are the number of nodes for XCOORD and YCOORD
The number of the connection pair is the same as the YCOORD set
The total series of the dataframe is
1 + 1 + 1 + countX *3 +countY * 3 + CountY
"""
import numpy as np
import pandas as pd
from math import *
import matplotlib.pyplot as plt


class DataPostProcess(object):
    def __init__(self, jobname, jobnumber, Stiffness, dx, InitialAngle):
        self.jobname = jobname
        self.jobnumber = jobnumber
        self.stiffness = Stiffness
        self.dx = dx
        self.InitialAngle = InitialAngle

    def Process(self):
        jobname = self.jobname
        jobnumber = self.jobnumber

        f = open('readmeStep.txt', 'r')
        step = f.readline()
        f.close()
        StepNumber = step

        f = open('Output_Connection.txt', 'r')
        connection = f.readline()
        f.close()
        ConnectionNumber = int(connection)

        Resultfilelist = []
        f = open('Resultfile_Label.txt', 'r')
        for line in f.readlines():
            Resultfilelist.append(line.strip("\n"))

        fpath = Resultfilelist[0] + ".txt"
        df = pd.read_csv(fpath)
        for i in range(1, len(Resultfilelist)):
            fpath_temp = Resultfilelist[i] + ".txt"
            df_temp = pd.read_csv(fpath_temp)
            df = pd.concat([df, df_temp], axis=1)

        resultname = "CreaseArray" + str(jobnumber) + str(StepNumber) + ".xlsx"
        df.to_excel(resultname, index=False)

        """
        Post process for the results from the excel file
        (1) The strain energy curves for plate and total;
        (2) gamma13 can be calculated from the crease coordinate;
        (3) Average Crease Angle for the all creases.
        """

        fpath = "CreaseArray" + str(jobnumber) + str(StepNumber) + ".xlsx"
        df = pd.read_excel(fpath)

        """
        (1) The strain energy curves for plate and total;
        """
        PlateALLSE = df.loc[:, "PlateALLSE"]
        TotalALLSE = df.loc[:, "TotalALLSE"]
        CreaseALLSE = TotalALLSE - PlateALLSE
        fig, ax = plt.subplots(2, 2, figsize=(8.5, 7))
        ax[0][0].set(xlabel='Time')
        ax[0][0].set(ylabel='ALLSE')
        ax[0][0].set_xlim(0, 1)

        CreaseData = pd.DataFrame()
        CreaseData['Time'] = df["time"]
        InitialEnergy = ConnectionNumber * self.stiffness * self.dx / 2.0 * 19 * (pi/2) ** 2
        Factor = self.stiffness * 100
        CreaseData['PlateAllSE'] = df["PlateALLSE"] / Factor
        CreaseData['TotalALLSE'] = (df["TotalALLSE"] + InitialEnergy) / Factor
        CreaseData['CreaseALLSE'] = (CreaseALLSE + InitialEnergy) / Factor
        ax[0][0].plot(CreaseData['Time'], CreaseData['PlateAllSE'], linewidth=2.0, label="PlateALLSE")
        ax[0][0].plot(CreaseData['Time'], CreaseData['TotalALLSE'], linewidth=2.0, label="TotalALLSE")
        ax[0][0].plot(CreaseData['Time'], CreaseData['CreaseALLSE'], linewidth=2.0, label="CreaseALLSE")
        ax[0][0].legend()

        """
        (2) gamma13 can be calculated from the crease coordinate;
        """
        DirectionData = pd.DataFrame()
        ArrayNumberY = int((ConnectionNumber + 2) / 3)
        for num in range(ArrayNumberY):
            OutPutXList = []
            for line in open("OutputX" + str(num + 1) + ".txt").readlines():
                OutPutXList.append(line.strip("\n"))
            NodeLabel = OutPutXList[0]
            NodeLabel = int(NodeLabel)
            vx1 = df["XCoordX." + str(NodeLabel)]
            vy1 = df["XCoordY." + str(NodeLabel)]
            vz1 = df["XCoordZ." + str(NodeLabel)]

            NodeLabel = OutPutXList[-1]
            NodeLabel = int(NodeLabel)
            vx2 = df["XCoordX." + str(NodeLabel)]
            vy2 = df["XCoordY." + str(NodeLabel)]
            vz2 = df["XCoordZ." + str(NodeLabel)]
            DirectionData["vx" + str(num + 1)] = vx2 - vx1
            DirectionData["vy" + str(num + 1)] = vy2 - vy1
            DirectionData["vz" + str(num + 1)] = vz2 - vz1

        for num in range(ArrayNumberY-1):
            x1 = DirectionData["vx"+str(num+1)].values
            x2 = DirectionData["vx"+str(num+2)].values
            y1 = DirectionData["vy"+str(num+1)].values
            y2 = DirectionData["vy"+str(num+2)].values
            z1 = DirectionData["vz"+str(num+1)].values
            z2 = DirectionData["vz"+str(num+2)].values

            temp1 = x1 * x2 + y1 * y2 + z1 * z2
            temp2 = np.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
            temp3 = np.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
            temp = np.arccos(temp1 / (temp2 * temp3))
            gamma13 = "gamma13" + str(num + 1) + "/pi"
            if num % 2 == 1:
                CreaseData[gamma13] = (-temp + pi) / pi
            else:
                CreaseData[gamma13] = (temp + pi) / pi

            ax[0][1].plot(CreaseData['Time'], CreaseData[gamma13], linewidth=2.0,
                          label=r'$\gamma_{13}-Vertex$'+str(num+1))

        ax[0][1].set(xlabel='Time')
        ax[0][1].set(ylabel=r'$\gamma_{13}$/$\pi$')
        ax[0][1].set_xlim(0, 1)
        ax[0][1].set_ylim(0, 2)
        ax[0][1].legend()

        """
        (3) Average Crease Angle for the all creases.
        """
        CUR1DataList = list(filter(lambda x: "Element ASSEMBLY" in x, df.columns.values.tolist()))
        ArrayNumberY = int((ConnectionNumber + 2) / 3)
        for num in range(ArrayNumberY):
            CreaseNumber = "Crease" + str(num + 1) + "/pi"
            countLY = len(open("OutputCrease" + str(num + 1) + ".txt").readlines())
            CUR1_temp = CUR1DataList[num * countLY: (num+1) * countLY]
            if num % 2 == 1:
                CreaseData[CreaseNumber] = (df[CUR1_temp].mean(axis=1))/pi
            else:
                CreaseData[CreaseNumber] = (df[CUR1_temp].mean(axis=1) + 2*pi)/pi

        ax[1][0].set(xlabel="time")
        ax[1][0].set(ylabel=r'$CreaseAngles/\pi$')
        # phi1 = CreaseData["Crease1/pi"]
        # phi2 = CreaseData["Crease2/pi"]
        # phi3 = CreaseData["Crease3/pi"]
        for num in range(ArrayNumberY):
            CreaseNumber = "Crease" + str(num + 1) + "/pi"
            ax[1][0].plot(CreaseData['Time'], CreaseData[CreaseNumber], linewidth=2.0, label="Crease" + str(num + 1) + "/pi")
        ax[1][0].set_xlim(0, 1)
        ax[1][0].set_ylim(0, 2)
        ax[1][0].legend()

        ax[1][1].set(xlabel='Time')
        ax[1][1].set(ylabel='ALLSE')
        ax[1][1].set_xlim(0, 1)
        ax[1][1].plot(CreaseData['Time'], CreaseData['TotalALLSE'], linewidth=2.0, label="TotalALLSE")
        ax[1][1].legend()


        figname = "MiuraArray-Job" + str(jobnumber) + str(StepNumber) +  "-Angle-Energy-Relationship"
        fig.suptitle(figname)
        fig.savefig(figname + '.jpg')
        plt.show()
        fig_name_xlsx = "CreaseData" + str(jobnumber) + str(StepNumber) + '.xlsx'
        CreaseData.to_excel(fig_name_xlsx, index=False)


        # """
        # (2) The X configuration for the specified frame;
        # Specified frame list Frames[] and the simulation process percent Prcess[]
        # """
        # CreaseXData = pd.DataFrame()
        # timelist = time.tolist()
        # timeFram = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        # Frames = []
        # for nu in range(len(timeFram)):
        #     answer = []
        #     for ni in timelist:
        #         answer.append(abs(timeFram[nu] - ni))
        #     Frames.append(answer.index(min(answer)))
        #
        # fig, ax = plt.subplots()
        # figname = "CreaseArray" + str(jobnumber) + str(StepNumber) + "XCoord"
        # fig.suptitle(figname)
        # plt.xlabel('Coordinate X')
        # plt.ylabel('Coordinate Z')
        # for i in range(len(Frames)):
        #     XListX = list(filter(lambda x: "XCoordX" in x, df.columns.values.tolist()))
        #     XListZ = list(filter(lambda x: "XCoordZ" in x, df.columns.values.tolist()))
        #     dfXListX = df[XListX].loc[Frames[i], :]
        #     dfXListZ = df[XListZ].loc[Frames[i], :]
        #     label = str(floor(timeFram[i] * 100)) + "%"
        #     CreaseXData['dfXListX' + str(i)] = dfXListX.tolist()
        #     CreaseXData['dfXListZ' + str(i)] = dfXListZ.tolist()
        #     ax.plot(dfXListX, dfXListZ, linewidth=2.0, label=label)
        # ax.legend()
        # fig.savefig(figname + '.jpg')
        # plt.show()
        # figname = "CreaseXData" + str(jobnumber) + str(StepNumber) + '.xlsx'
        # CreaseXData.to_excel(figname, index=False)


        # """
        # (3) The Y configuration for the specified frame;
        # """
        # CreaseYData = pd.DataFrame()
        # fig, ax = plt.subplots()
        # figname = "CreaseArray" + str(jobnumber)+ str(StepNumber) + "YCoord"
        # fig.suptitle(figname)
        # plt.xlabel('Coordinate Y')
        # plt.ylabel('Coordinate Z')
        # for i in range(len(Frames)):
        #     YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
        #     YListZ = list(filter(lambda x: "YCoordZ" in x, df.columns.values.tolist()))
        #     dfYListY = df[YListY].loc[Frames[i], :]
        #     dfYListZ = df[YListZ].loc[Frames[i], :]
        #     CreaseYData['dfYListY' + str(i)] = dfYListY.tolist()
        #     CreaseYData['dfYListZ' + str(i)] = dfYListZ.tolist()
        #     label = str(floor(timeFram[i] * 100)) + "%"
        #     ax.plot(dfYListY, dfYListZ, linewidth=2.0, label=label)
        # ax.legend()
        # fig.savefig(figname + '.jpg')
        # plt.show()
        #
        # """
        # (4) The Crease Spring Angle configuration for the specified frame;
        # """
        # fig, ax = plt.subplots()
        # figname = "CreaseArraye" + str(jobnumber) + str(StepNumber) + "CreaseMomentAngle"
        # fig.suptitle(figname)
        # plt.xlabel('Coordinate Y')
        # plt.ylabel('Crease Angle')
        # for i in range(len(Frames)):
        #     YListY = list(filter(lambda x: "YCoordY" in x, df.columns.values.tolist()))
        #     CMoment = list(filter(lambda x: "Element ASSEMBLY" in x, df.columns.values.tolist()))
        #     dfYListY = df[YListY].loc[Frames[i], :]
        #     dfCMoment = df[CMoment].loc[Frames[i], :]
        #     dfCangle = - dfCMoment / (self.stiffness * self.dx) + self.InitialAngle
        #     label = str(floor(timeFram[i] * 100)) + "%"
        #     CreaseYData['dfCangle' + str(i)] = dfCangle.tolist()
        #     ax.plot(dfYListY, dfCangle, linewidth=2.0, label=label)
        # ax.legend()
        # fig.savefig(figname + '.jpg')
        # plt.show()
        # figname = "CreaseYData" + str(jobnumber) + str(StepNumber) + '.xlsx'
        # CreaseYData.to_excel(figname, index=False)
        #
        # """
        # (5) The X configuration for the last frame
        # Specified frame list Frames[] and the simulation process percent Prcess[]
        # """
        # fig, ax = plt.subplots()
        # figname = "CreaseArray" + str(jobnumber) + str(StepNumber) +"XCoordLast"
        # fig.suptitle(figname)
        # plt.xlabel('Coordinate X')
        # plt.ylabel('Coordinate Z')
        # XListX = list(filter(lambda x: "XCoordX" in x, df.columns.values.tolist()))
        # XListZ = list(filter(lambda x: "XCoordZ" in x, df.columns.values.tolist()))
        # dfXListX = df[XListX].loc[df.index[-1], :]
        # dfXListZ = df[XListZ].loc[df.index[-1], :]
        # ax.plot(dfXListX, dfXListZ, linewidth=2.0)
        # fig.savefig(figname + '.jpg')
        # plt.show()

        # """
        # (6) The Y configuration parameter for the crease (average);
        # """
        # YListZ = list(filter(lambda x: "YCoordZ" in x, df.columns.values.tolist()))
        # dfYListZ = df[YListZ]
        # dfYListZNew = abs(dfYListZ)
        # NlistZpara = dfYListZNew.mean(axis=1)
        #
        # fig, ax = plt.subplots()
        # figname = "SingleCrease" + str(jobnumber) + "CreaseZ"
        # fig.suptitle(figname)
        # plt.xlabel('Time')
        # plt.ylabel('Coordinate Z Average')
        # ax.plot(time, NlistZpara, linewidth=2.0)
        # fig.savefig(figname + '.jpg')
        # plt.show()
        #
        # CreaseZ = pd.DataFrame()
        # CreaseZ['time'] = time.tolist()
        # CreaseZ['CreaseZ'] = NlistZpara.tolist()
        # figname = "CreaseZ" + str(jobnumber)+'.xlsx'
        # CreaseZ.to_excel(figname, index=False)


        # return AllSERatioResult
