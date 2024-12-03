#pragma once
#include "Common.h"




template<typename T>
class StencilWalk{
public:
    bool donorSearch(Point query_point, uint32 start_cell, uint32 &donor_cell);

    bool donorSearch(Point query_point, uint32 start_cell, uint32 &donor_cell, std::vector<uin32> &trace_path);
    StencilWalk(T &mesh_topology): mesh_topology_(mesh_topology);
private:
    // 搜索网格上的拓扑（网格单元包含的面
    const T &mesh_topology_;  

};







