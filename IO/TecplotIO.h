#pragma once
#include "Common.h"


class TecplotIO{
public:

    TecplotIO();
    void writeLineSegment(std::vector<Point> &points, bool independent_out=false, std::string file_name="");  // 线
    void writeFace();  // 面
    // 均以多面体方式写出
    void writeBlocks();  // 体
    std::ofstream *f_handle_;
private:
    TecplotIO(const TecplotIO &TecplotIO);  // override default copy constructor
    TecplotIO & operator = (const TecplotIO &TecplotIO);  // and assignment operator
}