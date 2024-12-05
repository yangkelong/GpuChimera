#include "TecplotIO.h"
#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */


// 即输出一系列点
// independent_out = false 时 将这部分写入的内容作为 .plt文件一部分
//                   true, 单独输出
void TecplotIO::writeLineSegment(std::vector<Point> &points, bool independent_out, std::string file_name){
    std::ofstream* f_handle = nullptr;
    std::string zone_name;
    if(independent_out){
        f_handle = new std::ofstream(file_name);
    }
    else
        f_handle = f_handle_;
    (*f_handle)<< "VARIABLES=\"X\",\"Y\",\"Z\" \n";
    (*f_handle)<< "ZONE T=\""<< zone_name << "\", I=" << points.size() << ", "<<"DATAPACKING=POINT\n";
    (*f_handle)<<std::setprecision(16)<<std::scientific;
    for(const auto& point: points){
        (*f_handle)<< point.x << "  " << point.y <<"  " <<  point.z <<"\n";
    }
    (*f_handle).flush();
    if(independent_out){
        f_handle->close();
        delete f_handle;
    }
}



void TecplotIO::writeBlocks(){
    // ?????????????? 索引应该是 0-based or 1-based, 面法向 面包含顶点 左手定则还是右手？？？, 左侧为
    std::string file_name;
    /* Call TECINI142 */
    INTEGER4 FileType = 0;   /* 0 for full file */
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT; Only PLT is currently
    INTEGER4 Debug = 0;
    INTEGER4 VIsDouble = 1;
    INTEGER4 res_code = 0;                      /* use to check return codes */
    // 初始化
    res_code = TECINI142((char *)"Pyramid",     /* Data Set Title */
                  (char *)"X Y Z iblank",       /* Variable List */
                //   (char *)"pyramid.plt", /* File Name */
                  file_name.data(), 
                  (char *)".",           /* Scratch Directory */
                  &FileFormat,
                  &(FileType),
                  &(Debug),
                  &(VIsDouble));
    INTEGER4 ZoneType = 7;   /* 7 for FEPolyhedron */
    for(const auto &block: blocks){
        INTEGER4 NumNodes = 5;  // vertex 数目
        INTEGER4 NumElems = 1;  // cell 数目
        INTEGER4 NumFaces = 5;  // facet 数目
        INTEGER4 ICellMax = 0;   /* Not Used, set to zero */
        INTEGER4 JCellMax = 0;   /* Not Used, set to zero */
        INTEGER4 KCellMax = 0;   /* Not Used, set to zero */
        double SolTime = 12.65;  /* solution time */
        INTEGER4 StrandID = 0;   /* static zone */
        INTEGER4 ParentZone = 0; /* no parent zone */
        INTEGER4 IsBlock = 1;    /* block format */
        INTEGER4 NFConns = 0;  /* not used for FEPolyhedron zones */
        INTEGER4 FNMode = 0;  /* not used for FEPolyhedron zones */
        INTEGER4 *PassiveVarArray = NULL;
        INTEGER4 *ValueLocationArray = {1, 1, 1, 0};  // 变量位置（格心 0, 格点 1, NULL=全部格点)  
        INTEGER4 *VarShareArray = NULL;
        INTEGER4 ShrConn = 0;
        /* The number of face nodes in the zone. This example creates
        * a zone with a single pyramidal cell. This cell has four
        * triangular faces and one rectangular face, yielding a total
        * of 16 face nodes.
        */
        INTEGER4 NumFaceNodes = 16;  // 全部 facet 包含顶点数目和
        INTEGER4 NumBConns = 0; /* No Boundary Connections */
        INTEGER4 NumBItems = 0; /* No Boundary Items */
        std::string &zone_name = block.zone_name;
        // 创建 zone
        res_code = TECZNE142((char *)"Polyhedral Zone (Octahedron)", &ZoneType, &NumNodes, &NumElems, &NumFaces,
                    &ICellMax, &JCellMax, &KCellMax, &SolTime, &StrandID, &ParentZone, &IsBlock,
                    &NFConns, &FNMode, &NumFaceNodes, &NumBConns, &NumBItems, PassiveVarArray,
                    ValueLocationArray, VarShareArray, &ShrConn);
        // 写入坐标
        double *X = new double[NumNodes];
        double *Y = new double[NumNodes];
        double *Z = new double[NumNodes];

        INTEGER4 DIsDouble = 1;  // 双精度浮点
        INTEGER4 DIsSingle = 0;  // 单精度浮点 or 整型？？？
        res_code = TECDAT142(&NumNodes, X, &DIsDouble);
        res_code = TECDAT142(&NumNodes, Y, &DIsDouble);
        res_code = TECDAT142(&NumNodes, Z, &DIsDouble);
        res_code = TECDAT142(&NumElems, blovk.iblanks, &DIsSingle);
        INTEGER4 *FaceNodeCounts = new INTEGER4[NumFaces];  // facet 包含顶点个数数组
        /* The first four faces are triangular, i.e. have three nodes.
        * The fifth face is rectangular, i.e. has four nodes. */
        FaceNodeCounts[0] = 3;
        FaceNodeCounts[1] = 3;
        FaceNodeCounts[2] = 3;
        FaceNodeCounts[3] = 3;
        FaceNodeCounts[4] = 4;
        INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];  // facet 包含顶点索引数组

    }   
    /* Face Nodes for Face 1 */
    FaceNodes[0] = 1;
    FaceNodes[1] = 2;
    FaceNodes[2] = 3;
    /* Face Nodes for Face 2 */
    FaceNodes[3] = 3;
    FaceNodes[4] = 2;
    FaceNodes[5] = 4;
    /* Face Nodes for Face 3 */
    FaceNodes[6] = 5;
    FaceNodes[7] = 2;
    FaceNodes[8] = 4;
    /* Face Nodes for Face 4 */
    FaceNodes[9] = 1;
    FaceNodes[10] = 2;
    FaceNodes[11] = 5;
    /* Face Nodes for Face 5 */
    FaceNodes[12] = 1;
    FaceNodes[13] = 5;
    FaceNodes[14] = 4;
    FaceNodes[15] = 3;
    
    INTEGER4 *FaceLeftElems = new INTEGER4[NumFaces];  // facet 左右网格单元索引
    FaceLeftElems[0] = 1;
    FaceLeftElems[1] = 1;
    FaceLeftElems[2] = 0;
    FaceLeftElems[3] = 0;
    FaceLeftElems[4] = 0;
    INTEGER4 *FaceRightElems = new INTEGER4[NumFaces];
    FaceRightElems[0] = 0;
    FaceRightElems[1] = 0;
    FaceRightElems[2] = 1;
    FaceRightElems[3] = 1;
    FaceRightElems[4] = 1;
    /* Write the face map (created above) using TECPOLYFACE142. */
    res_code = TECPOLYFACE142(&NumFaces, FaceNodeCounts, FaceNodes,
                              FaceLeftElems, FaceRightElems);
    delete FaceNodeCounts;
    delete FaceNodes;
    delete FaceLeftElems;
    delete FaceRightElems;
    // 关闭文件
    res_code = TECEND142();
    return 0;
}




