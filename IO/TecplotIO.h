#pragma once
#include "Common.h"


class TecplotIO{
public:

    TecplotIO();
    void readMesh();
    void constructIdNodes();
    vector<Point> & getCoordNodes();
    vector<NodeIdent> & getIdNodes();
    vector<NodeIdentMsh> & getIdNodesMsh();
    TecplotIO* setNbNode(unsigned);
    TecplotIO* setNbElMsh(unsigned);
    TecplotIO* setNbElm(unsigned);
    TecplotIO* setFileName(std::string);
    void writeLine();  // 线
    void writeFace();  // 面
    void writeBlock();  // 体

private:
    TecplotIO(const TecplotIO &TecplotIO);  // override default copy constructor
    TecplotIO & operator = (const TecplotIO &TecplotIO);  // and assignment operator
    unsigned nbNode;
    unsigned nbElMsh;
    unsigned nbElm;
    std::string fName;
    std::vector<Point> coordNodes;
    std::vector<NodeIdent> idNodes;
    std::vector<NodeIdentMsh> idNodesMsh;
}