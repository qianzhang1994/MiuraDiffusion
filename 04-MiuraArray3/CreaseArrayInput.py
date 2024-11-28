"""
Created on June 21 2022
Modified On November 26 2022
@author: Zhang Qian
"""

from math import *
import numpy as np


class CreaseArray(object):
    """
    Miura Array in the Y Direction.
    Loadregion, Creaseregion are not used in the model.
    """

    def __init__(self, jobname, Stiffness, dx, width, height, SectorAngle, InitialAngle, thickness, foldingState,
                 Loadregion, Creaseregion, ArrayNumberX, ArrayNumberY, CreaseNumber, CutNumber):
        self.jobname = jobname
        self.stiffness = Stiffness
        self.dx = dx
        self.width = width
        self.height = height
        self.SectorAngle = SectorAngle
        self.InitialAngle = InitialAngle
        self.thickness = thickness
        self.foldingstate = foldingState
        self.loadregion = Loadregion
        self.Creaseregion = Creaseregion
        self.ArrayNumberX = ArrayNumberX
        self.ArrayNumberY = ArrayNumberY
        self.CreaseNumber = CreaseNumber
        self.CutNumber = CutNumber

    def CreaseAngle(self, Angle1, Angle2, FoldingAngle):
        """"""
        CreaseAngle = acos(cos(Angle1)*cos(Angle2)+sin(Angle1)*sin(Angle2)*cos(FoldingAngle))
        return CreaseAngle

    def FoldAngle(self, Angle1, Angle2, Angle3):
        """"""
        FoldAngle = acos((cos(Angle3)-cos(Angle1)*cos(Angle2))/(sin(Angle1)*sin(Angle2)))
        return FoldAngle

    def _MiuraVertex(self, geometry, InitialAngle, ArrayNumberX, ArrayNumberY):
        """
        The vertex number of the plate starts from the left-up node,
        then increases count-clockwise, like
        1, 4
        2, 3
        height = distance between 1 and 2
        width = distance between 1 and 4
        The node 3 of left plate is selected as the origin point.
        the line 3/4 of left plate is the y-axis
        the 3/4 of left is connected by 2/1 of right plate
        :arg geometry [width, Height and sector angle] of the plate.
        :arg float InitialAngle: the initial angle of the two plate.
        :arg int
        """
        # ArrayNUmberX and ArrayNumberY is used for position
        width = geometry[0]
        height = geometry[1]
        SecAngle = geometry[2]
        CoordR = np.ones((4, 3))
        CreaseAngle = self.CreaseAngle(SecAngle, SecAngle, InitialAngle)
        CoordR[0] = np.array([0, 0, 0])
        CoordR[3] = np.array([width * sin(CreaseAngle/2), -width * cos(CreaseAngle/2), 0])

        FoldAngle = self.FoldAngle(SecAngle, CreaseAngle, SecAngle)
        Rotation = self.CreaseAngle(SecAngle, CreaseAngle/2, FoldAngle)
        temp = np.array([0, height * cos(Rotation), -height * sin(Rotation)])
        CoordR[1] = CoordR[0] + temp
        CoordR[2] = CoordR[3] + temp

        CoordArray = np.ones((4, 3))
        ArrayTempX = np.array([width * sin(CreaseAngle/2) * 2, 0, 0])
        ArrayTempY = np.array([0, 2 * height * cos(Rotation), 0])

        NumberX = (ArrayNumberX - 1) // 2
        NumberY = (ArrayNumberY - 1) // 2

        if (ArrayNumberX % 2 == 1) & (ArrayNumberY % 2 == 1):
            CoordArray = CoordR

        elif (ArrayNumberX % 2 == 0) & (ArrayNumberY % 2 == 1):
            CoordArray[0] = CoordR[3]
            CoordArray[1] = CoordR[2]
            CoordArray[2] = CoordR[1] + np.array([width * sin(CreaseAngle/2) * 2, 0, 0])
            CoordArray[3] = CoordR[0] + np.array([width * sin(CreaseAngle / 2) * 2, 0, 0])

        elif (ArrayNumberX % 2 == 1) & (ArrayNumberY % 2 == 0):
            CoordArray[0] = CoordR[1]
            CoordArray[1] = CoordR[0] + np.array([0, 2 * height * cos(Rotation), 0])
            CoordArray[2] = CoordR[3] + np.array([0, 2 * height * cos(Rotation), 0])
            CoordArray[3] = CoordR[2]

        elif (ArrayNumberX % 2 == 0) & (ArrayNumberY % 2 == 0):
            CoordArray[0] = CoordR[2]
            CoordArray[1] = CoordR[3] + np.array([0, 2 * height * cos(Rotation), 0])
            CoordArray[2] = CoordR[0] + np.array([width * sin(CreaseAngle/2) * 2, 2 * height * cos(Rotation), 0])
            CoordArray[3] = CoordR[1] + np.array([width * sin(CreaseAngle/2) * 2, 0, 0])

        CoordArray = CoordArray + NumberX * ArrayTempX + NumberY * ArrayTempY
        return CoordArray

    def _coordinate(self, geometry, InitialAngle, PlateNumber):
        """
        The vertex number of the plate starts from the left-up node,
        then increases count-clockwise, like
        1, 4
        2, 3
        The node 3 of left plate is selected as the origin point.
        the line 3/4 of left plate is the y-axis
        the 3/4 of left is connected by 2/1 of right plate
        :arg geometry [width, Height and sector angle] of the plate.
        :arg float InitialAngle: the initial angle of the two plate.
        :arg int PlateNumber: 0 for left plate and 1 for right plate
        """
        Rotation = (pi - InitialAngle) / 2
        T = np.array([[cos(Rotation), 0, sin(Rotation)],
                      [0, 1, 0],
                      [-sin(Rotation), 0, cos(Rotation)]])
        width = geometry[0]
        height = geometry[1]
        SecAngle = geometry[2]
        CoordR = np.ones((4, 3))

        CoordR[0] = np.array([0, height, 0])
        CoordR[1] = np.array([0, 0, 0])
        temp3 = np.array([width * sin(SecAngle), -width * cos(SecAngle), 0])
        CoordR[2] = np.dot(T, temp3)
        CoordR[3] = CoordR[2] + CoordR[0]

        if (PlateNumber % 2) == 1:
            Coord = CoordR
        else:
            CoordL = np.ones((4, 3))
            M = np.array([[-1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])
            CoordL[0] = np.dot(M, CoordR[3])
            CoordL[1] = np.dot(M, CoordR[2])
            CoordL[2] = CoordR[1]
            CoordL[3] = CoordR[0]
            Coord = CoordL

        # Translation Matrix in the Horizon Direction
        InitialWidth = width * sin(SecAngle) * cos(Rotation) * 2
        TM = np.array([InitialWidth, 0, 0])
        # The tension plane is defined as the coordinate plane
        InitialHeight = width * sin(SecAngle) * sin(Rotation)
        TZM = np.array([0, 0, InitialHeight])
        CoordArray = Coord + PlateNumber // 2 * TM + TZM

        return CoordArray

    def _Mesh(self, nw, nh, Coord, identifier):
        """
        Only support Quadrilateral mesh
        Calculate the node list and element list in global system
        Aspect_Ratio means the mesh size in the width and height direction.
        nw: width/dx/Aspect_Ratio
        nh: height/dx
        identifier: the plate number, int
        :return: NodeCoord Node list A(,4), (Node Number, x, y, z)
        :return: Element element lsit B(,5),
                (Element Number, Node1, Node2, Node 3, Node4 )
        :return: EdgeSet --Edge Node set for Connection design.
        Edge Number: 0 is for line 1-2, 1 is for line 2-3. 2 is for line 4-3
                3 is for line 1-4. The order of the line is considered.
        """
        NodeCoord = np.zeros(((nw + 1) * (nh + 1), 4))
        Element = np.zeros((nw * nh, 5), dtype=np.int32)
        EdgeSet = []
        # Node Information
        Number = 1
        for i in range(nw + 1):
            start = Coord[0] + (i / nw) * (Coord[3] - Coord[0])
            end = Coord[1] + (i / nw) * (Coord[2] - Coord[1])
            for j in range(nh + 1):
                temp = start + (j / nh) * (end - start)
                Global_Number = Number + identifier * (nw + 1) * (nh + 1)
                NodeCoord[Number - 1] = np.array([Global_Number, temp[0], temp[1], temp[2]])
                Number += 1
        # Element Information
        Number = 1
        for i in range(nw):
            for j in range(nh):
                Global_symbol = identifier * (nw + 1) * (nh + 1)
                Element[Number - 1] = np.array([Number + identifier * nw * nh,
                                                Global_symbol + i * (nh + 1) + j + 1,
                                                Global_symbol + i * (nh + 1) + j + 2,
                                                Global_symbol + (i + 1) * (nh + 1) + j + 2,
                                                Global_symbol + (i + 1) * (nh + 1) + j + 1])
                Number += 1
        # Edge Information for node set using the connection and loading
        for i in range(4):
            plateIncrease = identifier * (nw + 1) * (nh + 1)
            if (i % 2) == 0:
                EdgeIncrease = 1 if i == 0 else nw * (nh + 1) + 1
                EdgeNode = [itme for itme in range(plateIncrease + EdgeIncrease, plateIncrease + EdgeIncrease + nh + 1)]
                EdgeSet.append(EdgeNode)
            else:
                EdgeIncrease = nh + 1 if i == 1 else 1
                EdgeNode = [itme for itme in
                            range(plateIncrease + EdgeIncrease, plateIncrease + EdgeIncrease + (nh + 1) * nw + 1,
                                  nh + 1)]
                EdgeSet.append(EdgeNode)

        return NodeCoord, Element, EdgeSet

    def _Direction(self, Node_start, Node_end, Node_ref):
        x1 = Node_end - Node_start
        Direction1 = x1/np.linalg.norm(x1)
        x2 = Node_ref - Node_start
        Direction_temp = x2/np.linalg.norm(x2)
        Direction3 = np.cross(Direction1, Direction_temp)
        Direction2 = np.cross(Direction3, Direction1)
        Direction = np.concatenate((Direction1, Direction2))
        return Direction

    def Inpwrite(self):
        """
        Basic Information Input Module
        The plate geometry, [width, height, Sector Angle]
        the crease length is the height, and  Sector angle between line 1/4 to 1/2, default = 90
        InitialAngle: Define the initial rest angle of the one crease
        For the crease array, there are two kinds of plates, left and right.
        """
        width = self.width
        height = self.height
        SectorAngle = self.SectorAngle
        geometry = np.array([width, height, SectorAngle])
        InitialAngle = self.InitialAngle


        """
        Mesh Information
        dx is the size of mesh, nw, nh are the number of seed of
        edge with the width and height.
        NodeNumber and ElementNumber are the number of the nodes and 
        elements for each plate.
        dx > 10 * thickness for shell element
        Aspect Ratio is used to design the seed distribution in the width and height direction.
        If Aspect ratio is larger than 1.0, there are larger mesh size in the Height direction. 
        """
        dx = self.dx
        Aspect_Ratio = 1.0
        nw = floor(width / dx)
        nh = floor(height / dx / Aspect_Ratio)
        NodeNumber = (nw + 1) * (nh + 1)
        ElementNumber = nw * nh

        """
        Node List and Element List Modulus
        Assembly.
        Including NodeCoord and Element Information for all plate.
        The geometry points for each plate is defined as
        1  4
        2  3
        Plate Number: 0 is for left plate, while 1 is for right plate
        Edge Number: 0 is for line 1-2, 1 is for line 2-3. 2 is for line 4-3 
                    3 is for line 1-4. The order of the line is considered.
        Position of the edge node in the mesh, Plate Number, Edge Number, Node Number 
        such as Plate Number = 0, Edge Number = 1, Node Number = 5 means the nodes is located
        at the left plate, line 2-3, the fifth node. 

        """
        NodeCoord = []
        Element = []
        EdgeSet = []

        for j in range(self.ArrayNumberY):
            for i in range(self.ArrayNumberX):
                PlateNumber = j * self.ArrayNumberX + i
                Coord_temp = self._MiuraVertex(geometry, InitialAngle, i+1, j+1)
                NodeCoord_temp, Element_temp, Edge_temp = self._Mesh(nw, nh, Coord_temp, PlateNumber)
                if (i == 0) & (j == 0):
                    NodeCoord = NodeCoord_temp
                    Element = Element_temp
                    EdgeSet = Edge_temp
                else:
                    NodeCoord = np.concatenate((NodeCoord, NodeCoord_temp))
                    Element = np.concatenate((Element, Element_temp))
                    EdgeSet = EdgeSet + Edge_temp

        """
        Material Information of FEM Model
        The plate is elastic material
        Including Thickness, Elastic Module, Poisson ratio 
        unit [mm], [kg] [N] [rad]
        PET material except for the thickness
        ref: Embedded Actuation for Shape-Adaptive Origami, 2021, k=0.05
        B = Et^3/12(1-v^2)
        """
        thickness = self.thickness
        Elastic_Module = 3200.0
        Possion_ratio = 0.43
        Stiffness = self.stiffness
        CreaseStiffness = Stiffness * dx
        B = Elastic_Module * thickness ** 3 / 12 / (1 - Possion_ratio ** 2)
        ratio = Stiffness * self.width / B
        print("The scale parameter kL/B is ", ratio)

        Total_NodeNumber = self.ArrayNumberX * self.ArrayNumberY * NodeNumber
        Total_ElementNumber = self.ArrayNumberX * self.ArrayNumberY * ElementNumber




        """
        Connection Information
        Use the specified method instead of the query method
        For the two plate, the connection Information can be given 
        as (0,2) to (1,0) between the third edge of first plate and 
        the first edge of the second plate
        The local coordinate system is also defined through the node information
        There is only one connection edge
        ConectionEdge=[[2，4],[6，8],……]
        MVCrease , 1 is for mountain and -1 is for valley
        """
        ConnectionEdge = []
        ConnectionPlate = []
        ConnectionNumber = (self.ArrayNumberX - 1) * self.ArrayNumberY + self.ArrayNumberX * (self.ArrayNumberY - 1)
        MVCrease = []
        # X Connection
        for y_num in range(self.ArrayNumberY):
            for x_num in range(self.ArrayNumberX - 1):
                PlatePair = [x_num + y_num * self.ArrayNumberX, x_num + y_num * self.ArrayNumberX + 1]
                EdgePair = [2, 0]
                EdgePair_global = [4 * x + y for x, y in zip(PlatePair, EdgePair)]
                ConnectionEdge.append(EdgePair_global)
                ConnectionPlate.append([x_num + 1, y_num + 1, 2])
                if (x_num + y_num + 2) % 2 == 0:
                    MVCrease.append(1)
                else:
                    MVCrease.append(-1)

        for x_num in range(self.ArrayNumberX):
            for y_num in range(self.ArrayNumberY-1):
                PlatePair = [x_num + y_num * self.ArrayNumberX, x_num + (y_num + 1) * self.ArrayNumberX]
                EdgePair = [1, 3]
                EdgePair_global = [4 * x + y for x, y in zip(PlatePair, EdgePair)]
                ConnectionEdge.append(EdgePair_global)
                ConnectionPlate.append([x_num + 1, y_num + 1, 1])
                if (y_num + 1) % 2 == 0:
                    MVCrease.append(-1)
                else:
                    MVCrease.append(1)
        print("The Connection Pair (Edge Number) is ", ConnectionEdge)
        print("The Connection Pair (Plate Number) is ", ConnectionPlate)
        print("The Mountain and valley crease assignment is ", MVCrease)

        """
        Loading Modulus
        foldingState: 1 is unfolding and 0 is the folding process,default = 1
        if foldingstate = 0 and ArrayNumberY is odd (3, 5), U4 is positive (InitialAngle)
        if foldingstate = 0 and ArrayNumberY is even (2, 4), U4 is negative (-InitialAngle)
        if foldingstate = 1 and ArrayNumberY is odd (3, 5), U4 is negative (pi - InitialAngle)
        if foldingstate = 1 and ArrayNumberY is even (2, 4), U4 is positive (pi - InitialAngle)
        load boundary:
        Z   X   Z
            X
            X
        YZ  X   YZ
        NodeSetZ:  the U3 is constrained.
        NodeSetX:  the U1 is constrained.
        NodeSetYZ: the U2 and U3 are constrained.

        """
        if self.foldingstate == 0:
            if MVCrease[self.CreaseNumber - 1] == 1:
                Dis = pi - InitialAngle
            else:
                Dis = InitialAngle - pi
            Dis = 0.0
        else:
            if MVCrease[self.CreaseNumber - 1] == 1:
                Dis = InitialAngle - 1.0 * pi
            else:
                Dis = -InitialAngle + 1.0 * pi

        NodeSetZ = [1, (self.ArrayNumberX - 1) * (nw + 1) * (nh + 1) + nw * (nh + 1) + 1]
        NodeSetYZ1 = self.ArrayNumberX * (self.ArrayNumberY - 1) * (nw + 1) * (nh + 1) + (nh + 1)
        NodeSetYZ2 = self.ArrayNumberX * self.ArrayNumberY * (nw + 1) * (nh + 1)
        NodeSetYZ = [NodeSetYZ1, NodeSetYZ2]
        NodeSetX = []
        for i in range(self.ArrayNumberY):
            if i+1 not in self.CutNumber:
            # temp = self.ArrayNumberX * (self.ArrayNumberY - 1) * (nw + 1) * (nh + 1)
                temp = self.ArrayNumberX * i * (nw + 1) * (nh + 1)
                NodeSetX.append(temp + nw * (nh + 1) + 1)
                NodeSetX.append(temp + (nw + 1) * (nh + 1))
                NodeSetX.append(temp + (nw + 1) * (nh + 1) + 1)
                NodeSetX.append(temp + (nw + 1) * (nh + 1) + nh + 1)

        """
        Output Information for analysis
        Configuration Label
        OutputX : Coordinate of Creases.
        OutputLY OutputRY : Moment of Creases.
        """
        with open("Output_Connection.txt", "w") as f:
            f.write(str(ConnectionNumber) + "\n")

        for num in range(ConnectionNumber):
            EdgePair_global = ConnectionEdge[num]
            Pair1 = EdgePair_global[0]
            temp = EdgeSet[Pair1]
            OutputLY = temp[1:-1]
            with open("OutputCrease"+str(num + 1)+".txt", "w") as f:
                for i in range(len(OutputLY)):
                    f.write(str(OutputLY[i]) + "\n")

            with open("OutputX"+str(num + 1)+".txt", "w") as f:
                for i in range(len(temp)):
                    f.write(str(temp[i]) + "\n")


        """
        Inp file Module
        Label information: Width, Height, Sector Angle, Initial Angle, PlateNumber
        Define the name of the inp file as "CreaseArray"
        Static calculation in General/static solver. 
        """

        inp_file = open(self.jobname + ".inp", "w")
        inp_file.write("*Heading\n")
        inp_file.write("**Job Name and Model Name:CreaseArray\n")
        inp_file.write("*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        inp_file.write("**\n")

        inp_file.write("**PARTS\n")
        inp_file.write("*Part,name=Crease\n")
        inp_file.write("*Node\n")
        # for the node list
        for i in range(Total_NodeNumber):
            inp_file.write(str(i + 1) + "," + str(NodeCoord[i][1]) + ","
                           + str(NodeCoord[i][2]) + "," + str(NodeCoord[i][3]) + "\n")
        inp_file.write("**\n")

        inp_file.write("*Element,type=S4R\n")
        for i in range(Total_ElementNumber):
            inp_file.write(str(i + 1) + "," + str(Element[i][1]) + "," + str(Element[i][2]) + ","
                           + str(Element[i][3]) + "," + str(Element[i][4]) + "\n")
        inp_file.write("**\n")

        inp_file.write("*Nset, nset=SET-All, generate\n")
        inp_file.write(" 1," + str(Total_NodeNumber) + ",1\n")
        inp_file.write("*Elset, elset=SET-All, generate\n")
        inp_file.write("1," + str(Total_ElementNumber) + ",1\n")

        inp_file.write("** Section: \n")
        inp_file.write("*Shell Section, elset=SET-All, material=Self-define\n")
        inp_file.write(str(thickness) + ", 5\n")
        inp_file.write("*End Part\n")

        """
        Assembly Modulus, including the node set definition
        element set definition, and connection definition.
        """

        inp_file.write("** ASSEMBLY\n")
        inp_file.write("*Assembly, name=Assembly\n")
        inp_file.write("*Instance, name=CreaseArrayModel, part=Crease\n")
        inp_file.write("*End Instance\n")

        inp_file.write("*Nset, nset=Set-NodeZ, instance=CreaseArrayModel\n")
        inp_file.write(str(NodeSetZ[0]) + "," + str(NodeSetZ[1]) + "\n")

        inp_file.write("*Nset, nset=Set-NodeYZ, instance=CreaseArrayModel\n")
        inp_file.write(str(NodeSetYZ[0]) + "," + str(NodeSetYZ[1]) + "\n")

        inp_file.write("*Nset, nset=Set-NodeX, instance=CreaseArrayModel\n")
        for i in range(len(NodeSetX)):
            if i % 6 == 5:
                inp_file.write("\n")
            if i < len(NodeSetX) - 1:
                inp_file.write(str(NodeSetX[i]) + ",")
            else:
                inp_file.write(str(NodeSetX[i]) + "\n")

        inp_file.write("*Elset, elset=SET-Plate, instance=CreaseArrayModel, generate\n")
        inp_file.write("1," + str(Total_ElementNumber) + ",1\n")

        # there are many connection Pairs in this analysis
        # Each connection pair can have different definition.
        inp_file.write("*Element, type=CONN3D2\n")
        for num in range(ConnectionNumber):
            if num+1 not in self.CutNumber:
                EdgePair_global = ConnectionEdge[num]
                Pair1 = EdgePair_global[0]
                Pair2 = EdgePair_global[1]
                PairLen = len(EdgeSet[Pair1])
                # End Node are not coupled
                for i in range(1, PairLen-1):
                    if (i / PairLen < self.Creaseregion[0]) | (i / PairLen > self.Creaseregion[1]):
                        inp_file.write(str(i + (num * (PairLen - 2))) + ", CreaseArrayModel." + str(
                            EdgeSet[Pair1][i]) + ", CreaseArrayModel." + str(EdgeSet[Pair2][i]) + "\n")

        for num in range(ConnectionNumber):
            if num+1 not in self.CutNumber:
                EdgePair_global = ConnectionEdge[num]
                Pair1 = EdgePair_global[0]
                Pair2 = EdgePair_global[1]
                PairLen = len(EdgeSet[Pair1])
                inp_file.write("*Nset, nset=Set-Connection" + str(num + 1) + ", instance=CreaseArrayModel\n")
                for item in range(1, PairLen - 2):
                    if (item / PairLen < self.Creaseregion[0]) | (item / PairLen > self.Creaseregion[1]):
                        if item % 6 == 5:
                            inp_file.write("\n")
                        inp_file.write(str(EdgeSet[Pair1][item]) + "," + str(EdgeSet[Pair2][item]) + ",")
                inp_file.write(str(EdgeSet[Pair1][-2]) + "," + str(EdgeSet[Pair2][-2]) + "\n")
                inp_file.write("*Elset, elset=Set-Connection" + str(num + 1) + ", generate\n")
                inp_file.write(str(num * (PairLen - 2) + 1) + "," + str((num + 1) * (PairLen - 2)) + ",1 \n")

                # Orientation Calculation
                Plate_global = ConnectionPlate[num]
                Coord_Wire = self._MiuraVertex(geometry, InitialAngle, Plate_global[0], Plate_global[1])
                Node_start = Coord_Wire[Plate_global[2]]
                Node_end = Coord_Wire[(Pair1 + 1) % 4]
                Node_ref = Coord_Wire[(Pair1 + 2) % 4]
                Direction = self._Direction(Node_start, Node_end, Node_ref)

                inp_file.write("*Orientation, name=csys-Con" + str(num + 1) + "\n")
                inp_file.write(str(Direction[0])+","+str(Direction[1])+","+str(Direction[2])+","+str(Direction[3])+","
                               +str(Direction[4])+","+str(Direction[5])+"\n")
                inp_file.write("1,0.\n")
                inp_file.write("*Connector Section, elset=Set-Connection" + str(num + 1) + ", behavior=ConnProp-" +
                               str(num + 1) + "\n")
                inp_file.write("Hinge,\n")
                inp_file.write("csys-Con" + str(num + 1) + ",\n")
            # There are the same behavior for the hinge stiffness.
        inp_file.write("*End Assembly\n")
        for num in range(ConnectionNumber):
            if num+1 not in self.CutNumber:
                inp_file.write("*Connector Behavior, name=ConnProp-" + str(num + 1) + "\n")
                if MVCrease[num] == -1:
                    inp_file.write("*Connector Constitutive Reference\n")
                    inp_file.write(", , , 90., , \n")
                    inp_file.write("*Connector Elasticity, component=4\n")
                    if num < self.ArrayNumberY:
                        inp_file.write(str(CreaseStiffness) + ",\n")
                    else:
                        inp_file.write(str(CreaseStiffness/2.0) + ",\n")
                else:
                    inp_file.write("*Connector Constitutive Reference\n")
                    inp_file.write(", , , -90., , \n")
                    inp_file.write("*Connector Elasticity, component=4\n")
                    if num < self.ArrayNumberY:
                        inp_file.write(str(CreaseStiffness) + ",\n")
                    else:
                        inp_file.write(str(CreaseStiffness/2.0) + ",\n")

        """
        Standard Module
        Including Material Information, Step Information, 
        Boundary Information
        """
        inp_file.write("** MATERIALS\n")
        inp_file.write("*Material, name=Self-define\n")
        inp_file.write("*Elastic\n")
        inp_file.write(str(Elastic_Module) + "," + str(Possion_ratio) + "\n")

        inp_file.write("** STEP: Step-1\n")
        inp_file.write("*Step, name=Step-1, nlgeom=YES, inc=1500\n")
        inp_file.write("*Static, stabilize, factor = 1e-12, allsdtol = 0, continue=NO\n")
        # inp_file.write("*Static\n")
        inp_file.write("0.001, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-NodeZ, 3, 3 \n")
        inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-NodeX, 1, 1\n")
        inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-NodeYZ, 2, 2\n")
        inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        inp_file.write("*Boundary\n")
        inp_file.write("Set-NodeYZ, 3, 3, \n")
        # inp_file.write("** Name: BC-5 Type: Connector displacement\n")
        # inp_file.write("*Connector Motion\n")
        # inp_file.write("Set-Connection" + str(self.CreaseNumber) + ", 4, " + str(Dis) + "\n")


        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")

        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, number interval=10, time marks=YES\n")
        inp_file.write("** FIELD OUTPUT: F-Output-1\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-2\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection"+str(num+1)+", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-1\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-2\n")
        inp_file.write("*Output, history\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection"+str(num+1)+", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-3\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")

        # Setp-2
        inp_file.write("** STEP: Step-2\n")
        inp_file.write("*Step, name=Step-2, nlgeom=YES, inc=1500\n")
        inp_file.write("*Static, stabilize, factor = 1e-12, allsdtol = 0, continue=NO\n")
        # inp_file.write("*Static\n")
        inp_file.write("0.001, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        # inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeZ, 3, 3 \n")
        # inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeX, 1, 1\n")
        # inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 2, 2\n")
        # inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 3, 3, \n")
        inp_file.write("** Name: BC-5 Type: Connector displacement\n")
        inp_file.write("*Connector Motion\n")
        inp_file.write("Set-Connection" + str(self.CreaseNumber) + ", 4, " + str(Dis) + "\n")

        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")

        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, number interval=10, time marks=YES\n")
        inp_file.write("** FIELD OUTPUT: F-Output-3\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-4\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-4\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-5\n")
        inp_file.write("*Output, history\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-6\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")

        # Setp-3
        inp_file.write("** STEP: Step-3\n")
        inp_file.write("*Step, name=Step-3, nlgeom=YES, inc=1500\n")
        inp_file.write("*Static, stabilize, factor = 1e-12, allsdtol = 0, continue=NO\n")
        # inp_file.write("*Static\n")
        inp_file.write("0.001, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        # inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeZ, 3, 3 \n")
        # inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeX, 1, 1\n")
        # inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 2, 2\n")
        # inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 3, 3, \n")
        inp_file.write("** Name: BC-5 Type: Connector displacement\n")
        inp_file.write("*Connector Motion, op=NEW\n")
        # inp_file.write("Set-Connection" + str(self.CreaseNumber) + ", 4, " + str(Dis) + "\n")

        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")

        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, number interval=10, time marks=YES\n")
        inp_file.write("** FIELD OUTPUT: F-Output-3\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-4\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-4\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-5\n")
        inp_file.write("*Output, history\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-6\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")

        # Setp-4
        inp_file.write("** STEP: Step-4\n")
        inp_file.write("*Step, name=Step-4, nlgeom=YES, inc=1500\n")
        inp_file.write("*Static, stabilize, factor = 1e-12, allsdtol = 0, continue=NO\n")
        # inp_file.write("*Static\n")
        inp_file.write("0.001, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        # inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeZ, 3, 3 \n")
        # inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeX, 1, 1\n")
        # inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 2, 2\n")
        # inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 3, 3, \n")
        inp_file.write("** Name: BC-6 Type: Connector displacement\n")
        inp_file.write("*Connector Motion\n")
        Dis2 = 0.5
        inp_file.write("Set-Connection" + str(self.CreaseNumber) + ", 4, " + str(Dis2) + "\n")

        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")

        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, number interval=10, time marks=YES\n")
        inp_file.write("** FIELD OUTPUT: F-Output-3\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-4\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-4\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-5\n")
        inp_file.write("*Output, history\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-6\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")

        # Setp-5
        inp_file.write("** STEP: Step-5\n")
        inp_file.write("*Step, name=Step-5, nlgeom=YES, inc=1500\n")
        inp_file.write("*Static, stabilize, factor = 1e-12, allsdtol = 0, continue=NO\n")
        # inp_file.write("*Static\n")
        inp_file.write("0.001, 1., 1e-15, 0.01\n")

        inp_file.write("** BOUNDARY CONDITIONS\n")
        # inp_file.write("** Name: BC-1 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeZ, 3, 3 \n")
        # inp_file.write("** Name: BC-2 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeX, 1, 1\n")
        # inp_file.write("** Name: BC-3 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 2, 2\n")
        # inp_file.write("** Name: BC-4 Type: Displacement/Rotation\n")
        # inp_file.write("*Boundary\n")
        # inp_file.write("Set-NodeYZ, 3, 3, \n")
        inp_file.write("** Name: BC-6 Type: Connector displacement\n")
        inp_file.write("*Connector Motion, op=NEW\n")
        # inp_file.write("Set-Connection" + str(self.CreaseNumber) + ", 4, " + str(Dis) + "\n")

        inp_file.write("** CONTROLS\n")
        inp_file.write("*Controls, reset\n")
        inp_file.write("*Controls, parameters=time incrementation\n")
        inp_file.write("8, 10, , , , , , 50, , , \n")

        inp_file.write("** OUTPUT REQUESTS\n")
        inp_file.write("*Restart, write, number interval=10, time marks=YES\n")
        inp_file.write("** FIELD OUTPUT: F-Output-3\n")
        inp_file.write("*Output, field\n")
        inp_file.write("*Node Output\n")
        inp_file.write("CF, COORD, RF, U\n")
        inp_file.write("*Element Output, directions=YES\n")
        inp_file.write("LE, PE, PEEQ, PEMAG, S\n")
        inp_file.write("** FIELD OUTPUT: F-Output-4\n")
        inp_file.write("**\n")
        inp_file.write("*Output, field\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-4\n")
        inp_file.write("*Output, history, variable=PRESELECT\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-5\n")
        inp_file.write("*Output, history\n")
        for num in range(ConnectionNumber):
            if num + 1 not in self.CutNumber:
                inp_file.write("*Element Output, elset=Set-Connection" + str(num + 1) + ", directions=YES\n")
                inp_file.write("CTF, CU\n")
        inp_file.write("** HISTORY OUTPUT: H-Output-6\n")
        inp_file.write("*Output, history\n")
        inp_file.write("*Energy Output, elset=SET-Plate\n")
        inp_file.write("ALLSE, \n")
        inp_file.write("*End Step\n")

        inp_file.close()


