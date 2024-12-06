#pragma once
#include "Common.h"

struct TriFacet{
    uint32 id;  // 1-index
    uint32 parent_facet_id;  // 原始面元id
    uint32 twin_facet_id;  // 双胞胎三角面元id
    bool have_twin = false;  // 如果为四边形 facet 切分来, 具有一个双胞胎三角面元
    gdt::vec3i index;  // 3个顶点索引
    gdt::vec3f normal;  // 法向
    TriFacet(){}
    ~TriFacet() {}
};

class Block{
public:
    using Point = gdt::vec3f;
    int vertex_num, face_num, cell_num;
    int tet_num, pyramid_num, wedge_num, hex_num;
    int *cells_iblank;
    // 索引均保持 1-based
    Point *vertex_coord;
    // 面单元与顶点拓扑
    uint32 *face_ctn_vertex_xadj;
    uint32 *face_ctn_vertex_adjncy;
    gdt::vec2ui *face_side_cell;
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
    Point *cells_center;
    ~Block();

    /*  需要将四边形 facet 切分为 triangle facet
        并构建切分后的 拓扑
    */
    TriFacet *tri_facet_array;
    uint32 tri_facet_num;
    uint32 *cell_ctn_triface_xadj; 
    uint32 *cell_ctn_triface_adjncy;
    //
    void calCellCenter(); 
    void buildTriFacetArray();
    void getCellTriFacet(uint32 cell_id, std::set<uint32>&facet_set) const;
};