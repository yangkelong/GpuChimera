#pragma once
#include "../GridBlock.h"
#include "../Common.h"


class TecplotIO{
public:
    std::vector<Block*> blocks_ptr;
    TecplotIO()=default;
    // void writeLineSegment(std::vector<Point> &points, std::string file_name);  // 线
    void writeFace();  // 面
    void writeBlocks(const std::string &file_name);  // 体 均以多面体方式写出
};

