#include "GridBlock.h"

struct BoundaryCondition {
	const uint32 wall = 2;
};

struct FaceType {
	uint32 first_index;
	uint32 last_index;
	uint32 type;
	std::string name;
};

struct VolumeType {
	uint32 first_index;
	uint32 last_index;
	uint32 type;
	std::string name;
};

struct Cell {
	uint32 cellNo = 0;
	uint32 ctnVertexNum = 0;
	uint32 ctnVertexNo[8] = {};
	uint32 ctnFaceNum = 0;
	uint32 ctnFaceNo[6] = {};
	uint32 adjCellNum = 0;
	uint32 adjCellNo[6] = {};

	void setCellNo(uint32 cellNo) {
		this->cellNo = cellNo;
	}

	void addVertex(uint32 n, uint32 vertexArr[]) {
		if (cellNo == 0)
			return;
		const uint32 tmpN = this->ctnVertexNum;
		for (uint32 i = 0; i < n; ++i) {
			bool found = false;
			for (uint32 j = 0; j < tmpN; ++j) {
				if (vertexArr[i] == ctnVertexNo[j]) {
					found = true;
					break;
				}
			}
			// 不在当前数组中 则添加
			if (!found) {
				ctnVertexNo[ctnVertexNum] = vertexArr[i];
				++ctnVertexNum;
			}
		}
	}
	void addFace(uint32 faceNo) {
		if (cellNo == 0)
			return;
		this->ctnFaceNo[ctnFaceNum] = faceNo;
		++ctnFaceNum;
	}
	void addAdjCell(uint32 adjCellNo) {
		if (cellNo == 0)
			return;
		this->adjCellNo[adjCellNum] = adjCellNo;
		++adjCellNum;
	}
};

    Block::Block(const char *file_name){
        readMesh(file_name);
    }

    Block::~Block(){
        delete []vertex_coord;
        delete []face_side_cell;
        delete []cell_ctn_vertex_xadj;
        delete []cell_ctn_vertex_adjncy;
        delete []cell_adj_cell_xadj;
        delete []cell_adj_cell_adjncy;
        delete []cell_ctn_face_xadj;
        delete []cell_ctn_face_adjncy;
        delete []face_ctn_vertex_xadj;
        delete []face_ctn_vertex_adjncy;
        delete []bc_faces;
        delete []interior_faces;
        delete []cells_iblank;
    }

// 网格读取成功返回0
int Block::readMesh(const char *meshfilename){
    uint32 *faceCtnVertexAdjncy;
    uint32 n0, n1, n2, n3, c0, c1;
    uint32 faceCtnVertex[4];
    std::vector<FaceType> FaceTypeSet;
    std::vector<VolumeType> VolumeZoneSet;
    std::vector<Cell> cellSet;
    orderDict zoneID2index;
    unsigned int index, i, VertexNum;
    char ch, buff[1000];
    char bctype[50], facename[50];
    int DIM;
    unsigned int first_index, last_index;
    unsigned int iVertex, iFace;
    uint32 zone_id, type, shapetype, cell_type;
    int count = 0;     //	用于index = 45时，统计网格面类型的计数器
    int kkk, OKOK = 0; //	用于处理index ＝0的情况
    int offset;
    std::cout << "\n>>Reading gambit mesh file, please wait......" <<std::string(meshfilename)<< std::endl;
    std::ifstream file_handle(meshfilename);
    if (!file_handle) //	文件打开出错，返回-1
    {
        std::cout << "\n\aERROR: Unable to open " << meshfilename << " for reading!\n\a";
        std::cin.get();
        return -1;
    }
    while (!file_handle.eof())
    {
        // 读直到遇到字符 '('
        do
        {
            file_handle.get(ch);
            if (file_handle.eof())
                break;
        } while (ch != '(');
        if (file_handle.eof())
            break;
        OKOK = 0;
        // 读取 字符'(' 后的十六进制数据转化为十进制 并存入 index
        file_handle >> std::dec >> index;
        // std::cout << index << "start!!!" << std::endl;
        switch (index)
        {
        // 0 注释  (0 "comment text")
        case 0:
            kkk = 1;
            do
            {
                file_handle.get(ch);
                if (ch == '(')
                    kkk++;
                if (ch == ')')
                    kkk--;
                if (ch == '"')
                {   
                    file_handle.getline(buff, 500, '"');
                    std::cout << buff << std::endl;
                    std::string str(buff);
                    const std::string to_find_tet = "Tet cells";
                    const std::string to_find_pyramid = "Pyramid cells";
                    const std::string to_find_wedge = "Wedge cells";
                    const std::string to_find_hex = "Hex cells";
                    const std::string to_find_bc = "Boundary Faces";
                    const std::string to_find_interior = "Interior Faces";
                    std::regex pattern(R"((\d+))"); // 匹配一个或多个数字
                    std::smatch matches;

                    std::regex_search(str, matches, pattern);
                    const int number = std::stoi(matches[1].str());

                    if (str.find(to_find_tet) != std::string::npos){
                        tet_num = number;
                    }
                    else if(str.find(to_find_pyramid) != std::string::npos){
                        pyramid_num = number;
                    }
                    else if(str.find(to_find_wedge) != std::string::npos){
                        wedge_num = number;
                    }
                    else if(str.find(to_find_hex) != std::string::npos){
                        hex_num = number;
                    }
                    else if(str.find(to_find_bc) != std::string::npos){
                        bc_face_num = number;
                        bc_faces = new uint32[bc_face_num];
                    }
                    else if(str.find(to_find_interior) != std::string::npos){
                        interior_face_num = number;
                        interior_faces = new uint32[interior_face_num];
                    }
                }
            } while (kkk != 0);
            OKOK = 1;
            break;
        // 2 维度 (2 ND)
        case 2:
            file_handle >> DIM;
            if (DIM != 3)
            {
                std::cout << "current not support 2D mesh";
                return 1;
            }
            break;
        // 10 顶点  (10 (zone-id first-index last-index type ND)( x1 y1 z1 ...
        case 10:
        case 2010:
        case 3010:
            do{
                file_handle.get(ch);
            } while (ch != '(');

            file_handle >> std::hex >> zone_id >> first_index >> last_index >> type;  //  ND;
            do{
                file_handle.get(ch);
            } while (ch != ')');
            // 总的网格顶点信息
            if (zone_id == 0)
            {
                std::cout << "Total grid vertexs:" << std::setw(8) << last_index << std::endl;
                vertex_num = last_index;
                vertex_coord = new float64[3*(vertex_num+1)];
            }
            // 每个zone的网格点信息
            if (zone_id != 0)
            {
                do{
                    file_handle.get(ch);
                } while (ch != '(');
                for (iVertex = first_index; iVertex <= last_index; iVertex++){
                    // 读网格点坐标
                    file_handle >> std::dec >> vertex_coord[iVertex * 3] >> vertex_coord[iVertex * 3 + 1] >> vertex_coord[iVertex * 3 + 2];
                }
                do
                {
                    file_handle.get(ch);
                } while (ch != ')');
            }
            break;
        //	12 网格 (12 (zone-id first-index last-index type element-type))
        case 12:
        case 2012:
        case 3012:
            do
            {
                file_handle.get(ch);
            } while (ch != '(');

            file_handle >> std::hex >> zone_id;
            //	网格单元的总的信息
            if (zone_id == 0)
            {
                file_handle >> std::hex >> first_index >> last_index >> shapetype;
                do
                {
                    file_handle.get(ch);
                } while (ch != ')');
                std::cout << "Total grid cells:" << std::setw(8) << last_index << std::endl;
                cell_num = last_index;
                cellSet.resize(cell_num + 1);
                cells_iblank = new int[cell_num + 1];
                assert(cell_num == tet_num + pyramid_num + wedge_num + hex_num);
                cell_ctn_vertex_xadj = new uint32[cell_num + 2];
                cell_ctn_vertex_adjncy = new uint32[tet_num * 4 + pyramid_num * 5 + wedge_num * 6 + hex_num * 8 + 1];
                cell_ctn_face_adjncy = new uint32[tet_num * 4 + pyramid_num * 5 + wedge_num * 6 + hex_num * 8 + 1];
                cell_adj_cell_adjncy = new uint32[tet_num * 4 + pyramid_num * 5 + wedge_num * 6 + hex_num * 8 + 1];
                cell_adj_cell_xadj = new uint32[cell_num + 2];
                cell_ctn_face_xadj = new uint32[cell_num + 2];
                cell_ctn_face_xadj[1] = 1;
                cell_adj_cell_xadj[1] = 1;
                cell_ctn_vertex_xadj[1] = 1;
                for (uint32 i = 1; i <= cell_num; ++i)
                    cellSet[i].setCellNo(i);
            }
            //	每个zone的网格单元信息
            if (zone_id != 0)
            {
                VolumeType volumeType;
                // FileRead >> std::hex >> first_index >> last_index >> type >> shapetype;
                file_handle >> std::hex >> volumeType.first_index >> volumeType.last_index >> volumeType.type >> shapetype;
                //	储存该体计算域的信息
                VolumeZoneSet.push_back(volumeType);
                zoneID2index[zone_id] = (int)VolumeZoneSet.size() - 1;
                do
                {
                    file_handle.get(ch);
                } while (ch != ')');
                if (shapetype == 0) //	混合网格
                {
                    do
                    {
                        file_handle.get(ch);
                    } while (ch != '(');
                    for (i = first_index; i <= last_index; i++)
                        file_handle >> std::hex >> cell_type; // 网格单元的形状，本程序暂时没用到
                    do
                    {
                        file_handle.get(ch);
                    } while (ch != ')');
                }
            }
            break;
        //	13 面   (13 (zone-id first-index last-index bc-type face-type))
        case 13:
        case 2013:
        case 3013:
            do
            {
                file_handle.get(ch);
            } while (ch != '(');

            file_handle >> std::hex >> zone_id;
            if (zone_id == 0) //	总的网格面信息
            {
                file_handle >> std::hex >> first_index >> last_index >> shapetype;
                do
                {
                    file_handle.get(ch);
                } while (ch != ')');
                std::cout << "Total grid faces:" << std::setw(8) << last_index << std::endl;
                face_num = last_index;
                face_side_cell = new uint32[2*(face_num + 1)];
                face_ctn_vertex_xadj = new uint32[face_num + 2];
                faceCtnVertexAdjncy = new uint32[face_num*4 + 1];
                face_ctn_vertex_xadj[1] = 1;

            }
            else if (zone_id != 0) //	每类型网格面的信息
            {
                FaceType faceType;
                // FileRead >> std::hex >> first_index >> last_index >> type >> shapetype;
                file_handle >> std::hex >> faceType.first_index >> faceType.last_index >> faceType.type >> shapetype;
                do
                {
                    file_handle.get(ch);
                } while (ch != ')');
                //		设置网格面的类型
                FaceTypeSet.push_back(faceType);
                zoneID2index[zone_id] = (unsigned int)FaceTypeSet.size() - 1; //	建立zone_id和类型的映射关系，后面读取边界

                do
                {
                    file_handle.get(ch);
                } while (ch != '(');

                VertexNum = shapetype;
                // if (shapetype == 5)
                //	IsPolyhedral = true;
                for (iFace = faceType.first_index; iFace <= faceType.last_index; iFace++)
                {
                    if (shapetype == 0 || shapetype == 5)
                        file_handle >> std::hex >> VertexNum; //	顶点个数
                    // 读取面单元数据   x n0 n1 ... nf c0 c1
                    // Note: 此处不作边界面处的处理 默认边界面应该是满足 c1 == 0， 在完成拓扑数组后针对 所有边界面检查一下是否满足 c1 == 0
                    if (VertexNum == 3){
                        file_handle >> std::hex >> faceCtnVertex[0] >> faceCtnVertex[1] >> faceCtnVertex[2] >> c0 >> c1;
                        uint32 tmpIndex = face_ctn_vertex_xadj[iFace];
                        face_ctn_vertex_xadj[iFace + 1] = tmpIndex + 3;
                        faceCtnVertexAdjncy[tmpIndex] = faceCtnVertex[2];
                        faceCtnVertexAdjncy[tmpIndex + 1] = faceCtnVertex[1];
                        faceCtnVertexAdjncy[tmpIndex + 2] = faceCtnVertex[0];
                        face_side_cell[iFace * 2] = c0;
                        face_side_cell[iFace * 2 + 1] = c1;
                        cellSet[c0].addVertex(3, faceCtnVertex);
                        cellSet[c0].addFace(iFace);
                        cellSet[c0].addAdjCell(c1);
                        cellSet[c1].addVertex(3, faceCtnVertex);
                        cellSet[c1].addFace(iFace);
                        cellSet[c1].addAdjCell(c0);
                    }
                    else if (VertexNum == 4){
                        file_handle >> std::hex >> faceCtnVertex[0] >> faceCtnVertex[1] >> faceCtnVertex[2] >> faceCtnVertex[3] >> c0 >> c1;
                        // faceCtnVertexXadj[iFace + 1] = faceCtnVertexXadj[iFace] + 4;
                        uint32 tmpIndex = face_ctn_vertex_xadj[iFace];
                        face_ctn_vertex_xadj[iFace + 1] = tmpIndex + 4;
                        faceCtnVertexAdjncy[tmpIndex] = faceCtnVertex[3];
                        faceCtnVertexAdjncy[tmpIndex + 1] = faceCtnVertex[2];
                        faceCtnVertexAdjncy[tmpIndex + 2] = faceCtnVertex[1];
                        faceCtnVertexAdjncy[tmpIndex + 3] = faceCtnVertex[0];
                        face_side_cell[iFace * 2] = c0;
                        face_side_cell[iFace * 2 + 1] = c1;
                        cellSet[c0].addVertex(4, faceCtnVertex);
                        cellSet[c0].addFace(iFace);
                        cellSet[c0].addAdjCell(c1);
                        cellSet[c1].addVertex(4, faceCtnVertex);
                        cellSet[c1].addFace(iFace);
                        cellSet[c1].addAdjCell(c0);
                    }
                    else
                        return -2;
                    // std::cout << VertexNum<<" "<< c0 << "  " << c1 << std::endl;
                }
                do
                {
                    file_handle.get(ch);
                } while (ch != ')');
            }
            // std::cout << "13 面" << "done!!!" << std::endl;
            break;
        // (39 (zone-id zone-type zone-name domain-id)(
        //     (condition1.value1)
        case 39:
        case 2039:
        case 3039:
            do
            {
                file_handle.get(ch);
            } while (ch != '(');

            file_handle >> zone_id >> bctype;
            file_handle.ignore(5, ' ');
            for (i = 0; i < 50; i++)
                facename[i] = '\0';
            offset = 0;
            file_handle.get(ch);
            while (ch != ' ' && ch != ')')
            {
                facename[offset] = ch;
                offset++;
                file_handle.get(ch);
                if (ch == ' ')
                {
                    file_handle.getline(buff, 100, ')');
                    break;
                }
            }
            file_handle.putback(')');
            do
            {
                file_handle.get(ch);
            } while (ch != ')');
            if (strcmp(bctype, "fluid\0") == 0)
            {
                index = zoneID2index[zone_id];
                VolumeZoneSet[index].name = facename;
            }
            else
            { //	面的信息
                index = zoneID2index[zone_id];
                FaceTypeSet[index].name = facename;
            }
            do
            {
                file_handle.get(ch);
            } while (ch != '(');
            do
            {
                file_handle.get(ch);
            } while (ch != ')');
            break;
        // (45 (zone-id zone-type zone-name domain-id)(
        //     (condition1.value1)
        case 45: //	Gambit和TGrid和pointwise(名字正序）
        case 2045:
        case 3045:
            do
            {
                file_handle.get(ch);
            } while (ch != '(');

            file_handle >> zone_id >> bctype;
            file_handle.ignore(5, ' ');
            file_handle.getline(facename, 100, ')');
            file_handle.putback(')');
            do
            {
                file_handle.get(ch);
            } while (ch != ')');
            // 流体域
            if (strcmp(bctype, "fluid\0") == 0)
            {
                index = zoneID2index[zone_id];
                VolumeZoneSet[index].name = facename;
            }
            else
            { //	面的信息
                index = zoneID2index[zone_id];
                FaceTypeSet[index].name = facename;
            }
            do
            {
                file_handle.get(ch);
            } while (ch != '(');
            do
            {
                file_handle.get(ch);
            } while (ch != ')');
            break;
        default: //
            kkk = 1;
            do
            {
                file_handle.get(ch);
                if (ch == '(')
                    kkk++;
                if (ch == ')')
                    kkk--;
            } while (kkk != 0);
            OKOK = 1;
            break;
        }
        if (OKOK == 1)
            continue;
        do
        {
            file_handle.get(ch);
            if (file_handle.eof())
                break;
        } while (ch != ')');
        if (file_handle.eof())
            break;
    }
    file_handle.close();
    // cellCtnVertex, cellCtnFace, cellAdjCell
    for (uint32 cellNo = 1; cellNo <= cell_num; ++cellNo)
    {
        Cell &cell = cellSet[cellNo];
        // build cellCtnVertex
        uint32 tmpIndex = cell_ctn_vertex_xadj[cellNo];
        cell_ctn_vertex_xadj[cellNo + 1] = tmpIndex + cell.ctnVertexNum;
        for (uint32 i = 0; i < cell.ctnVertexNum; ++i)
            cell_ctn_vertex_adjncy[tmpIndex + i] = cell.ctnVertexNo[i];
        // build cellCtnFace
        tmpIndex = cell_ctn_face_xadj[cellNo];
        cell_ctn_face_xadj[cellNo + 1] = tmpIndex + cell.ctnFaceNum;
        for (uint32 i = 0; i < cell.ctnFaceNum; ++i)
            cell_ctn_face_adjncy[tmpIndex + i] = cell.ctnFaceNo[i];
        // build cellAdjCell
        tmpIndex = cell_adj_cell_xadj[cellNo];
        cell_adj_cell_xadj[cellNo + 1] = tmpIndex + cell.adjCellNum;
        for (uint32 i = 0; i < cell.adjCellNum; ++i)
            cell_adj_cell_adjncy[tmpIndex + i] = cell.adjCellNo[i];
    }
    // 边界面, 内部面
    const uint32 interior = 2;
    const uint32 wall = 3;
    const uint32 pressure_inlet = 4;
    const uint32 pressure_outlet = 5;
    const uint32 symmetry = 7;
    const uint32 pressure_far_field = 9;
    // (BCsInfo["wall"].shape[0], BCsInfo["symmetry"].shape[0], BCsInfo["farField"].shape[0], BCsInfo["pressureInlet"].shape[0], BCsInfo["pressureOutlet"].shape[0])
    mylist wallFaces;
    mylist symmetryFaces;
    mylist farFieldFaces;
    mylist pressureInletFaces;
    mylist pressureOutletFaces;
    mylist interiorFaces;
    uint32 bc_faces_inf[5];

    for (FaceType &faceType : FaceTypeSet)
    {
        switch (faceType.type)
        {
        case interior:
            for (uint32 index = faceType.first_index; index <= faceType.last_index; ++index)
                interiorFaces.push_back(index);
            break;
        case wall:
            bc_faces_inf[0] += faceType.last_index - faceType.first_index + 1;
            for (uint32 index = faceType.first_index; index <= faceType.last_index; ++index)
                wallFaces.push_back(index);
            break;
        case pressure_inlet:
            bc_faces_inf[3] += faceType.last_index - faceType.first_index + 1;
            for (uint32 index = faceType.first_index; index <= faceType.last_index; ++index)
                pressureInletFaces.push_back(index);
            break;
        case pressure_outlet:
            bc_faces_inf[4] += faceType.last_index - faceType.first_index + 1;
            for (uint32 index = faceType.first_index; index <= faceType.last_index; ++index)
                pressureOutletFaces.push_back(index);
            break;
        case symmetry:
            bc_faces_inf[1] += faceType.last_index - faceType.first_index + 1;
            for (uint32 index = faceType.first_index; index <= faceType.last_index; ++index)
                symmetryFaces.push_back(index);
            break;
        case pressure_far_field:
            bc_faces_inf[2] += faceType.last_index - faceType.first_index + 1;
            for (uint32 index = faceType.first_index; index <= faceType.last_index; ++index)
                farFieldFaces.push_back(index);
            break;
        }
    }
    index = 0;
    for (auto faceNo : wallFaces)
    {
        bc_faces[index] = faceNo;
        wall_faces.push_back(faceNo);
        ++index;
    }
    for (auto faceNo : symmetryFaces)
    {
        bc_faces[index] = faceNo;
        ++index;
    }
    for (auto faceNo : farFieldFaces)
    {
        bc_faces[index] = faceNo;
        ++index;
    }
    for (auto faceNo : pressureInletFaces)
    {
        bc_faces[index] = faceNo;
        ++index;
    }
    for (auto faceNo : pressureOutletFaces)
    {
        bc_faces[index] = faceNo;
        ++index;
    }
    const uint32 bcFaceSize = index;
    index = 0;
    for (auto faceNo : interiorFaces){
        interior_faces[index] = faceNo;
        ++index;
    }
    // 完成拓扑数组后针对 所有边界面检查一下是否满足 c1 == 0
    for (index = 0; index < bcFaceSize; ++index)
    {
        uint32 faceNo = bc_faces[index];
        if (face_side_cell[faceNo * 2 + 1] != 0)
            return -3;
    }
    // 
    const int var = face_ctn_vertex_xadj[face_num + 1];
    face_ctn_vertex_adjncy = new uint32 [var];
    for (int i = 0; i < var; ++i){
        face_ctn_vertex_adjncy[i] =faceCtnVertexAdjncy[i];
        // std::cout<<"?????"<<faceCtnVertexAdjncy[i]<<std::endl;
    }
    delete []faceCtnVertexAdjncy;
    return 0;
}
