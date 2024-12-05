#include "TecplotIO.h"
#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */


// 即输出一系列点
// void TecplotIO::writeLineSegment(std::vector<Point> &points, std::string file_name){
//     std::ofstream* f_handle = nullptr;
//     std::string zone_name;
//     if(independent_out){
//         f_handle = new std::ofstream(file_name);
//     }
//     else
//         f_handle = f_handle_;
//     (*f_handle)<< "VARIABLES=\"X\",\"Y\",\"Z\" \n";
//     (*f_handle)<< "ZONE T=\""<< zone_name << "\", I=" << points.size() << ", "<<"DATAPACKING=POINT\n";
//     (*f_handle)<<std::setprecision(16)<<std::scientific;
//     for(const auto& point: points){
//         (*f_handle)<< point.x << "  " << point.y <<"  " <<  point.z <<"\n";
//     }
//     (*f_handle).flush();
//     if(independent_out){
//         f_handle->close();
//         delete f_handle;
//     }
// }


void TecplotIO::writeBlocks(const std::string &file_name){
    // 索引 1-based, 面包含顶点 右手定则 指向R侧 c1
    // 对于 <=01111111 11111111 11111111 11111111(2^31-1) 的uint32 其编码与 int32相同, 此处直接将 uint32 数组写入
    /* Call TECINI142 */
    INTEGER4 FileType = 0;   /* 0 for full file */
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT; Only PLT is currently
    INTEGER4 Debug = 0;
    INTEGER4 VIsDouble = 1;
    INTEGER4 res_code = 0;                      /* use to check return codes */
    INTEGER4 ZoneType = 7;   /* 7 for FEPolyhedron */
    int block_index = 0;
    // 初始化
    res_code = TECINI142((char *)"XXX",     /* Data Set Title */
               (char *)"X Y Z iblank",       /* Variable List */
            //    (char *)"X Y Z",       /* Variable List */
               file_name.data(), 
               (char *)".",           /* Scratch Directory */
               &FileFormat,
               &(FileType),
               &(Debug),
               &(VIsDouble));
    for(const auto &block_ptr: blocks_ptr){
       INTEGER4 NumNodes = block_ptr->vertex_num;  // vertex 数目
       INTEGER4 NumElems = block_ptr->cell_num;  // cell 数目
       INTEGER4 NumFaces = block_ptr->face_num;  // facet 数目
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
       INTEGER4 ValueLocationArray[4] = {1, 1, 1, 0};  // 变量位置（格心 0, 格点 1, NULL=全部格点)  
    //    INTEGER4 *ValueLocationArray = NULL;  // 变量位置（格心 0, 格点 1, NULL=全部格点)  
       
       INTEGER4 *VarShareArray = NULL;
       INTEGER4 ShrConn = 0;
       /* The number of face nodes in the zone. This example creates
       * a zone with a single pyramidal cell. This cell has four
       * triangular faces and one rectangular face, yielding a total
       * of 16 face nodes.
       */
       INTEGER4 NumFaceNodes = (INTEGER4)block_ptr->face_ctn_vertex_xadj[NumFaces+1]-1;  // 全部 facet 包含顶点数目和
       std::cout << "NumNodes: " << NumNodes << std::endl;
       std::cout << "NumElems: " << NumElems << std::endl;
       std::cout << "NumFaces: " << NumFaces << std::endl;
       INTEGER4 NumBConns = 0; /* No Boundary Connections */
       INTEGER4 NumBItems = 0; /* No Boundary Items */
       std::string zone_name("zone_" + std::to_string(block_index));
       // 创建 zone
       res_code = TECZNE142((char *)"Polyhedral Zone (Octahedron)", &ZoneType, &NumNodes, &NumElems, &NumFaces,
                   &ICellMax, &JCellMax, &KCellMax, &SolTime, &StrandID, &ParentZone, &IsBlock,
                   &NFConns, &FNMode, &NumFaceNodes, &NumBConns, &NumBItems, PassiveVarArray,
                   ValueLocationArray, VarShareArray, &ShrConn);
       // 写入坐标
       double *X = new double[NumNodes];
       double *Y = new double[NumNodes];
       double *Z = new double[NumNodes];
       for(int i=1; i<=NumNodes; ++i){
           X[i-1] = block_ptr->vertex_coord[3*i];
           Y[i-1] = block_ptr->vertex_coord[3*i+1];
           Z[i-1] = block_ptr->vertex_coord[3*i+2];
       }
       INTEGER4 DIsDouble = 1;  // 双精度浮点
       INTEGER4 DIsSingle = 0;  // 单精度浮点 or 整型？？？
       for(int i=1; i<=NumElems; ++i){
        block_ptr->cells_iblank[i] = i%2;
       }
       res_code = TECDAT142(&NumNodes, X, &DIsDouble);
       res_code = TECDAT142(&NumNodes, Y, &DIsDouble);
       res_code = TECDAT142(&NumNodes, Z, &DIsDouble);
       res_code = TECDAT142(&NumElems, block_ptr->cells_iblank +1, &DIsSingle);
       INTEGER4 *FaceNodeCounts = new INTEGER4[NumFaces];  // facet 包含顶点个数数组
       for(int i=1; i<=NumFaces; ++i){
           FaceNodeCounts[i-1] = block_ptr->face_ctn_vertex_xadj[i+1]-block_ptr->face_ctn_vertex_xadj[i];
       }
       INTEGER4 *FaceNodes = reinterpret_cast<INTEGER4*>(block_ptr->face_ctn_vertex_adjncy+1);  // facet 包含顶点索引数组
       for(int i=1; i<=NumFaces; ++i){
           //for(int s=block_ptr->face_ctn_vertex_xadj[i]; s<block_ptr->face_ctn_vertex_xadj[i+1]; ++s)
           //    std::cout<<FaceNodes[s]<<" ";
           //std::cout<<std::endl;
       }
       INTEGER4 *FaceLeftElems = new INTEGER4[NumFaces];  // facet 左右网格单元索引
       INTEGER4 *FaceRightElems = new INTEGER4[NumFaces];
       for(int i=1; i<=NumFaces; ++i){
            FaceLeftElems[i-1] = block_ptr->face_side_cell[2*i];
            FaceRightElems[i-1] = block_ptr->face_side_cell[2*i+1];
            // std::cout << FaceLeftElems[i - 1] << " || " << FaceRightElems[i - 1] << std::endl;
       }
       /* Write the face map (created above) using TECPOLYFACE142. */
       res_code = TECPOLYFACE142(&NumFaces,
                       FaceNodeCounts,  /* The face node counts array */
                       FaceNodes,       /* The face nodes array */
                       FaceLeftElems,   /* The left elements array */
                       FaceRightElems); /* The right elements array */
       delete []FaceNodeCounts;
       delete []FaceLeftElems;
       delete []FaceRightElems;
       block_index += 1;
    }   
    // 关闭文件
    res_code = TECEND142();
    return;
}




