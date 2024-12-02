#define PROGRAM_ONCE
#include "Common.h"


class Block{
public:
    int vertex_num, face_num, cell_num;
    int tet_num, pyramid_num, wedge_num, hex_num;
    // 索引均保持 1-based
    float64 *vertex_coord;
    // 面单元与顶点拓扑
    uint32 *face_ctn_vertex_xadj;
    uint32 *face_ctn_vertex_adjncy;
    uint32 *face_side_cell;
    // 网格单元与顶点拓扑
    uint32 *cell_ctn_vertex_xadj;
    uint32 *cell_ctn_vertex_adjncy;
    // 网格单元与面拓扑
	uint32 *cell_ctn_face_xadj; 
    uint32 *cell_ctn_face_adjncy;
    // 网格单元与网格单元拓扑
    uint32 *cell_adj_cell_xadj;
    uint32 *cell_adj_cell_adjncy; 
    // 边界面（嵌套网格中只需要区别 wall no-wall）
    std::vector<uint32> wall_faces;
    uint32 *bc_faces;
    uint32 *interior_faces;
    int bc_face_num, interior_face_num;
    Block(const char *);
    int readMesh(const char *);
    ~Block();

};