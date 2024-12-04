from numba import cuda
from math import sqrt, ceil
import os
from DubheGpuSolver.setting import profile
from DubheGpuSolver.setting.profile import OutputDir, TecplotDataPrecision, THREADSPERBLOCK, INT32, lib, \
    MPI_isRoot, MPI_Root, UINT32, pltFluidDataSuffix, pltTopoDataSuffix, nVar, PRECISION, Crt_LENGTH
from DubheGpuSolver.init import clock
from DubheGpuSolver.Solver.Solver import conVarBar2fluidVar
from enum import IntEnum
from DubheGpuSolver.PostProcess.tecio import *
from DubheGpuSolver.Communication.MyMPI import MyMPI, MPI_Size
from DubheGpuSolver.PostProcess.dataSerialization import loadData, saveData
from DubheGpuSolver.MeshSystem.Geometry import coordScaling


class FieldDataType(IntEnum):
    Float = 1
    Double = 2
    Int32 = 3
    Int16 = 4
    Byte = 5


def _ctring(s):
    return ctypes.c_char_p(s.encode())


def _list_to_int_array(l):
    t = (ctypes.c_int * len(l))(*l)
    return t


# Note: 当前数据在写入plt 文件前均为双精度( MPI通信数据也全都是双精度), 只有在写入plt时 转成了float32类型写入
class TecplotIO:
    # 输出类初始化时将网格信息保存, 避免后续每次调用时反复处理拓扑信息, 当前只保存合并后的zone
    def __init__(self, vertices_coord, amountInfo, faces_sideCellsNo, faceCtnVertex, BCsInfo, debug=True, isRoot=False):
        # 流场变量个数
        self.nVar = 9
        if isRoot:
            if not profile.pltGenerateOnline:
                data = [vertices_coord, amountInfo, faces_sideCellsNo, faceCtnVertex, BCsInfo]
                saveData(data, file=OutputDir + profile.gridName + pltTopoDataSuffix)
                self.amountInfo = amountInfo
                return
            # 坐标为无量纲量需要转为有量纲
            n = vertices_coord.shape[0] - 1
            scaleFactor = Crt_LENGTH
            vertices_coord_device = cuda.to_device(vertices_coord)
            coordScaling[ceil(n / THREADSPERBLOCK), THREADSPERBLOCK](n, vertices_coord_device, scaleFactor)
            cuda.synchronize()
            self.vertices_coord = vertices_coord_device.copy_to_host()
            self.faces_sideCellsNo = faces_sideCellsNo
            self.coord = np.ascontiguousarray(self.vertices_coord[1:].T.flatten())
            # 后续这一块每一个zone都单独需要
            # num_nodes, num_faces, num_elements
            self.amountInfo = amountInfo
            faceCtnVertexXadj, faceCtnVertexAdjncy = faceCtnVertex[0], faceCtnVertex[1]
            self.total_num_face_nodes = faceCtnVertexAdjncy.shape[0]-1
            self.face_nodes = faceCtnVertexAdjncy[1:]
            num_faces = faceCtnVertexXadj.shape[0]-2
            self.face_node_counts = np.zeros((num_faces, ), dtype=INT32)
            lib.c_buildFaceNodeCounts(ctypes.c_int32(num_faces), faceCtnVertexXadj.ctypes, self.face_node_counts.ctypes)
            faces_sideCellsNo = faces_sideCellsNo.T
            self.left_elems = np.ascontiguousarray(faces_sideCellsNo[0, 1:])
            self.right_elems = np.ascontiguousarray(faces_sideCellsNo[1, 1:])
            self.valueLocation = [1] * 3 + [0] * self.nVar
            if TecplotDataPrecision == np.float32:
                self.VIsDouble = 0
            else:
                self.VIsDouble = 1
            self.debug = debug
            # 0 .plt    1 .szplt
            self.fileFormat = 0
            self.boundaryInfo = dict()
            # 初始化构建壁面等边界面 拓扑 ; 当前只有壁面 作为单独zone 保存
            wall_faces = BCsInfo["wall"]
            self.initialBoundaryTopo(wall_faces, faceCtnVertex)

    def initialBoundaryTopo(self, bcFaces, faceCtnVertex):
        uint32 = ctypes.c_uint32
        faceCtnVertexXadj, faceCtnVertexAdjncy = faceCtnVertex[0], faceCtnVertex[1]
        bcFacesSize = bcFaces.shape[0]
        # 该zone 包含的顶点个数 边个数; 结构 非结构网格以四边形估计最大值
        # maxNum_nodes = 4*bcFacesSize
        # maxNum_edges = 4*bcFacesSize
        # 多面体网格估计最大值; 如果面不相邻，则顶点个数最大为faceCtnVertexAdjncy.shape[0]; 边个数同样最大为该值
        maxNum_nodes = faceCtnVertexAdjncy.shape[0]
        maxNum_edges = faceCtnVertexAdjncy.shape[0]
        # 以边界面包含的顶点数量最大值 maxNum_nodes 初始化, 建立拓扑后需调整大小
        vertexNoArr_Global = np.zeros((maxNum_nodes, ), dtype=UINT32)
        face_nodes = np.zeros((2*maxNum_edges, ), dtype=UINT32)
        left_elems = np.zeros((maxNum_edges, ), dtype=UINT32)
        right_elems = np.zeros((maxNum_edges, ), dtype=UINT32)
        num_nodes = uint32(0)
        num_faces = uint32(0)
        num_elements = uint32(0)
        lib.c_buildBoundaryTopo(uint32(bcFacesSize), bcFaces.ctypes, faceCtnVertexXadj.ctypes, faceCtnVertexAdjncy.ctypes,
                                vertexNoArr_Global.ctypes, ctypes.byref(num_nodes), ctypes.byref(num_faces),
                                ctypes.byref(num_elements), face_nodes.ctypes, left_elems.ctypes, right_elems.ctypes)
        num_nodes = num_nodes.value
        num_faces = num_faces.value
        num_elements = num_elements.value
        vertexNoArr_Global = vertexNoArr_Global[:num_nodes]
        coord = np.ascontiguousarray(self.vertices_coord[vertexNoArr_Global].T.flatten())
        face_nodes = face_nodes[:2*num_faces]
        left_elems = left_elems[:num_faces]
        right_elems = right_elems[:num_faces]
        total_num_face_nodes = face_nodes.shape[0]
        # 边界流场数据从边界面左侧即内侧单元获取, lCellsNo-1 得到0起始的索引数组
        lCellsNo = self.faces_sideCellsNo[bcFaces][:, 0]
        indexArr = lCellsNo - 1
        # offset = scale * cellsAmount
        indexArr = np.concatenate([indexArr + scale*self.amountInfo[2] for scale in range(self.nVar)])
        self.boundaryInfo["wall"] = [num_nodes, num_faces, num_elements, total_num_face_nodes, coord, face_nodes, left_elems,
                                     right_elems, indexArr]
        pass

    def write_szplt(self, dataset):
        self._writerOpen(dataset)
        zones = dataset['zones']
        for zone_name in zones.keys():
            data = zones[zone_name]
            # xyz格点  vars格心
            nodes_xyz = data['xyz']
            zoneTitle = _ctring(data['zoneTitle'])
            #
            zoneType = data['zoneType']
            vars_data = data['vars']
            n = vars_data.shape[0]
            numberOfNodes = data['numberOfNodes']
            numberOfCells = data['numberOfCells']
            variableDataTypes = _list_to_int_array(data['variableDataTypes'])
            shareVarFromZone = None
            valueLocation = _list_to_int_array(data['valueLocation'])
            passiveVarList = None
            shareConnectivityFromZone = None
            numFaceConnections = 0
            faceNeighborMode = 2
            zoneIndex = ctypes.c_int(0)
            ref_zoneIndex = ctypes.byref(zoneIndex)
            tecio.tecZoneCreateFE(self.filehandle, zoneTitle, zoneType, numberOfNodes, numberOfCells,
                                     variableDataTypes, shareVarFromZone, valueLocation, passiveVarList,
                                     shareConnectivityFromZone, numFaceConnections, faceNeighborMode, ref_zoneIndex)
            solutionTime = 0
            strandID = 0
            self._tecZoneSetUnsteadyOptions(zoneIndex, solutionTime, strandID)
            # 变量数据n + 坐标数据3
            for k in range(n + 3):
                variableIndex = k + 1
                if k < 3:
                    var_data = nodes_xyz[k, :]
                else:
                    var_data = vars_data[k - 3, :]
                fieldDataType = variableDataTypes[variableIndex - 1]
                if fieldDataType == FieldDataType.Float:
                    tecio.tecZoneVarWriteFloatValues(self.filehandle, zoneIndex, variableIndex, 0, var_data.size,
                                                        var_data.ctypes)
                elif fieldDataType == FieldDataType.Double:
                    tecio.tecZoneVarWriteDoubleValues(self.filehandle, zoneIndex, variableIndex, 0,
                                                         var_data.size, var_data.ctypes)
                elif fieldDataType == FieldDataType.Int32:
                    tecio.tecZoneVarWriteInt32Values(self.filehandle, zoneIndex, variableIndex, 0, var_data.size,
                                                        var_data.ctypes)
                elif fieldDataType == FieldDataType.Int16:
                    tecio.tecZoneVarWriteInt16Values(self.filehandle, zoneIndex, variableIndex, 0, var_data.size,
                                                        var_data.ctypes)
                elif fieldDataType == FieldDataType.Byte:
                    tecio.tecZoneVarWriteUint8Values(self.filehandle, zoneIndex, variableIndex, 0, var_data.size,
                                                        var_data.ctypes)
                else:
                    raise Exception('FieldDataType Error:not defined data type')
            nodeMap = data['nodeMap']
            numValues = nodeMap.size
            # nodemap 中的节点编号 std::vector<int64_t> nodeMap(numValues);如果是int64 调用write64; .itemsize  元素字节数
            if nodeMap.itemsize == 8:
                tecio.tecZoneNodeMapWrite64(self.filehandle, zoneIndex, 0, 1, numValues, nodeMap.ctypes)
            else:
                tecio.tecZoneNodeMapWrite32(self.filehandle, zoneIndex, 0, 1, numValues, nodeMap.ctypes)
        tecio.tecFileWriterClose(self.ref_filehandle)

    def output_contour_bin(self, cells_W, nPointDomain, fname_contour='contour.plt', gridType=1):
        start = clock()
        # 传递无量纲的守恒量至根进程, 根进程将这些数据序列化到文件
        if not profile.pltGenerateOnline:
            if MPI_Size > 1:
                # 组装发送数组
                sendArr = np.ascontiguousarray(cells_W[1:nPointDomain+1].copy_to_host().T.flatten())
                itemSize = nVar
                # 根进程接收全部数据
                recvArr = MyMPI.gatherCellDataForTecplot(sendArr, itemSize, MPI_Root)
                if MPI_isRoot:
                    recvArr = recvArr.reshape((nVar, self.amountInfo[2])).T
                    cells_W = np.concatenate([np.zeros((1, nVar)), recvArr])
                # cells_W = np.zeros(self.amountInfo[2])
            else:
                cells_W = cells_W.copy_to_host()
            if MPI_isRoot:
                data = [profile.Crt_VELOCITY, profile.Crt_PRESSURE, profile.Ref_PRESSURE, profile.Ref_DynamicPressure, cells_W]
                saveData(data, OutputDir+fname_contour[:-4]+pltFluidDataSuffix)
            return
        data_device = cuda.to_device(np.zeros((self.nVar*nPointDomain), dtype=PRECISION))
        # 无量纲的守恒量转化为有量纲的需要输出的变量
        conVarBar2fluidVar[ceil(nPointDomain / THREADSPERBLOCK), THREADSPERBLOCK](nPointDomain, cells_W, data_device)
        cuda.synchronize()
        data = data_device.copy_to_host()
        if MPI_Size > 1:
            # 子进程传递数据至根进程, 根进程完成tecplot文件写入
            # 组装发送数组
            sendArr = data
            itemSize = self.nVar
            # 根进程接收全部数据
            recvArr = MyMPI.gatherCellDataForTecplot(sendArr, itemSize, MPI_Root)
            if MPI_isRoot:
                data = recvArr
            else:
                return
        # "dataType": 变量默认格式  float
        dataset = {"fileName": OutputDir+fname_contour,
                   "dataSetTitle": "FE zone",
                   "listOfVariables": ["X", "Y", "Z", "Density", "Vx", "Vy", "Vz", "V",  "M", "Pressure",
                                       "Temperature", "Cp"],
                   "fileType": 0,
                   "dataType": 1,
                   "zones": {1: {"zoneTitle": "zone1", "zoneType": 4,
                                 'vars': data, "variableDataTypes": [2] * 11, "valueLocation": self.valueLocation,
                                 }, }
                   }
        if fname_contour.endswith(".szplt"):
            # .szplt 导出时需要基于cell的拓扑
            raise Exception("current not support")
            nodeMap = data['nodeMap']
            cells_node = None
            # .cas 文件 cells_ctnVerticesNo,cells_node为一个列表, 元素为集合
            if gridType == 0:
                cells_node = cells_node[1:]
                cells_node = [list(item) for item in cells_node]
                cells_node = np.array(cells_node).flatten()
            # .cas 文件 cells_ctnVerticesNo 为 (cellCtnVertexXadj, cellCtnVertexAdjncy)
            elif gridType == 1:
                cells_node = cells_node[1][1:]
            self.write_szplt(dataset)
        else:
            self.write_plt(dataset)
        print("write done! spend time: %f" % (clock() - start))

    def write_plt(self, dataset):
        FName = dataset['fileName']
        Title = dataset['dataSetTitle']
        Variables = dataset['listOfVariables']
        Debug = self.debug
        open_file(FName, Title, Variables, self.VIsDouble, Debug)
        fluidData = None
        # 流场数据
        zones = dataset['zones']
        for zone_name in zones.keys():
            data = zones[zone_name]
            num_nodes, num_faces, num_elements = self.amountInfo[0], self.amountInfo[1], self.amountInfo[2]
            total_num_face_nodes = self.total_num_face_nodes
            face_nodes = self.face_nodes
            left_elems = self.left_elems
            right_elems = self.right_elems
            face_node_counts = self.face_node_counts
            # 创建 zone
            create_poly_zone(str(zone_name), ZONETYPE_FEPOLYHEDRON, num_nodes, num_elements, num_faces,
                             total_num_face_nodes, value_locations=data["valueLocation"])
            # 写入 坐标
            zone_write_arr(self.coord, self.VIsDouble)
            # zone_write_values(self.coord)
            # 写入 流场数据
            fluidData = data['vars']
            zone_write_arr(data['vars'], self.VIsDouble)
            # facet 数目, facet 包含顶点个数数组, facet 包含顶点索引数组, facet 左右网格单元索引
            tecpolyface(num_faces, face_node_counts, face_nodes, left_elems, right_elems)
        # 边界数据, 每个边界作为单独的zone输出
        # 需要更改, 如果流场数据按zone分开输出，当前 zones只有一个zone，整个流场数据只包含在一个zone中
        self.writeBoundary(fluidData)
        close_file()

    def writeBoundary(self, fluidData):
        for zone_name in self.boundaryInfo.keys():
            num_nodes, num_faces, num_elements, total_num_face_nodes, coord, face_nodes, left_elems, \
                right_elems, indexArr = self.boundaryInfo[zone_name]
            create_poly_zone(str(zone_name), ZONETYPE_FEPOLYGON, num_nodes, num_elements, num_faces,
                             total_num_face_nodes, value_locations=self.valueLocation)
            # 写入 坐标
            zone_write_arr(coord, self.VIsDouble)
            data = fluidData[indexArr]
            # 写入 流场数据
            zone_write_arr(data, self.VIsDouble)
            tecpolyface(num_faces, None, face_nodes, left_elems, right_elems)
        pass

    def _writerOpen(self, dataset):
        self.filehandle = ctypes.pointer(ctypes.c_int(0))
        self.ref_filehandle = ctypes.byref(self.filehandle)
        filename = _ctring(dataset['fileName'])
        dataSetTitle = _ctring(dataset['dataSetTitle'])
        listOfVariables = _ctring(",".join(dataset['listOfVariables']))
        if dataset['fileName'].endswith('.szplt'):
            fileFormat = 1
        else:
            raise Exception('file format error, not a .szplt format')
        fileType = dataset['fileType']
        dataType = dataset['dataType']
        gridFileHandle = None
        tecio.tecFileWriterOpen(filename, dataSetTitle, listOfVariables, fileFormat, fileType, dataType,
                                   gridFileHandle, self.ref_filehandle)
        if self.debug:
            outputDebugInfo = 1
            tecio.tecFileSetDiagnosticsLevel(self.filehandle, outputDebugInfo)

    def _tecZoneSetUnsteadyOptions(self, zone_n, solutionTime=0, StrandID=0):
        if solutionTime != 0 or StrandID != 0:
            solutionTime = ctypes.c_double(solutionTime)
            tecio.tecZoneSetUnsteadyOptions(self.filehandle, zone_n, solutionTime, StrandID)


# 云服务器上保存单元守恒量等数据
# 本地生成.plt文件
def generatePltOnLocal(gridName):
    vertices_coord, amountInfo, faces_sideCellsNo, faceCtnVertex, BCsInfo = loadData(file=OutputDir + gridName + pltTopoDataSuffix)
    tecplot = TecplotIO(vertices_coord, amountInfo, faces_sideCellsNo, faceCtnVertex, BCsInfo, isRoot=True)
    # 读取将该文件夹下所有的 pltFluidDataSuffix 后缀的文件, 生成对应的 .plt文件
    # Crt_VELOCITY, Crt_PRESSURE, Ref_PRESSURE, Ref_DynamicPressure
    nPoint = amountInfo[2]
    for filename in os.listdir(OutputDir):
        if filename.endswith(pltFluidDataSuffix):
            # loadData
            Crt_VELOCITY, Crt_PRESSURE, Ref_PRESSURE, Ref_DynamicPressure, cells_W = loadData(OutputDir+filename)
            np.save(gridName+".npy", cells_W)
            profile.Crt_VELOCITY = Crt_VELOCITY
            profile.Crt_PRESSURE = Crt_PRESSURE
            profile.Ref_PRESSURE = Ref_PRESSURE
            profile.Ref_DynamicPressure = Ref_DynamicPressure
            #
            filename = filename[:-len(pltFluidDataSuffix)] + ".plt"
            tecplot.output_contour_bin(cuda.to_device(cells_W), nPoint, filename)
    pass


if __name__ == '__main__':
    # profile.OutputDir = ""
    profile.pltGenerateOnline = True
    generatePltOnLocal("ball300")
    # for filename in os.listdir(OutputDir):
    #     if filename.endswith(pltFluidDataSuffix):
    #         # loadData
    #         filename = filename[:-len(pltFluidDataSuffix)] + ".plt"
    #         print(filename)
    pass
