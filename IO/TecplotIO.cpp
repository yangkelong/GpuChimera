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
        INTEGER4 NumNodes = 5;   /* number of unique nodes */
        INTEGER4 NumElems = 1;   /* number of elements */
        INTEGER4 NumFaces = 5;   /* number of unique faces */
        INTEGER4 ICellMax = 0;   /* Not Used, set to zero */
        INTEGER4 JCellMax = 0;   /* Not Used, set to zero */
        INTEGER4 KCellMax = 0;   /* Not Used, set to zero */
        double SolTime = 12.65;  /* solution time */
        INTEGER4 StrandID = 0;   /* static zone */
        INTEGER4 ParentZone = 0; /* no parent zone */
        INTEGER4 IsBlock = 1;    /* block format */
        INTEGER4 NFConns = 0;    /* not used for FEPolyhedron
                                * zones
                                */
        INTEGER4 FNMode = 0;     /* not used for FEPolyhedron
                                * zones
                                */
        INTEGER4 *PassiveVarArray = NULL;
        /*The location of each variable in the data set. ValueLocation(I) indicates the location of variable I
        for this zone. 0=cell-centered, 1=node-centered. Pass null to indicate that all variables are node-centered.*/
        INTEGER4 *ValueLocationArray = NULL;  
        INTEGER4 *VarShareArray = NULL;
        INTEGER4 ShrConn = 0;
        /* The number of face nodes in the zone. This example creates
        * a zone with a single pyramidal cell. This cell has four
        * triangular faces and one rectangular face, yielding a total
        * of 16 face nodes.
        */
        INTEGER4 NumFaceNodes = 16;
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
        INTEGER4 DIsDouble = 1; /* One for double precision */
        res_code = TECDAT142(&NumNodes, X, &DIsDouble);
        res_code = TECDAT142(&NumNodes, Y, &DIsDouble);
        res_code = TECDAT142(&NumNodes, Z, &DIsDouble);
        INTEGER4 *FaceNodeCounts = new INTEGER4[NumFaces];
        /* The first four faces are triangular, i.e. have three nodes.
        * The fifth face is rectangular, i.e. has four nodes. */
        FaceNodeCounts[0] = 3;
        FaceNodeCounts[1] = 3;
        FaceNodeCounts[2] = 3;
        FaceNodeCounts[3] = 3;
        FaceNodeCounts[4] = 4;
        INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];
        // 流场数据

    }   


    // 边界数据, 每个边界作为单独的zone输出
}




