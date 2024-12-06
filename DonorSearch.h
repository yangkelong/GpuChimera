#pragma once
#include "Common.h"
#include "GridBlock.h"


class StencilWalk {
public:
    using Point = Block::Point;
    uint32 max_step = 10;
    bool donorSearch(const Point &query_point, uint32 start_cell, uint32& donor_cell);
    bool donorSearch(const Point &query_point, uint32 start_cell, uint32& donor_cell, std::vector<Point> &trace_path);
    StencilWalk(Block& mesh_) : mesh(mesh_) {};
private:
    const Block& mesh;

};