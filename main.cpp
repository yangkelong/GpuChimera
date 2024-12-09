#include "GridBlock.h"
#include "IO/TecplotIO.h"
#include "DonorSearch.h"

// #include "TecplotWriter.h"


int main(int argc, char *argv[]) {
    Block mesh(argv[1]);
    TecplotIO tec_out;
    Block::Point query_point(2.1, 2.9, 2.5);
    uint32 start_cell = 3;
    uint32 donor_cell = 0;
    std::vector<Block::Point> trace_path;
    //trace_path.emplace_back(0,0,0);
    //trace_path.emplace_back(1,1,0);
    //trace_path.emplace_back(1,0,0);
    //trace_path.emplace_back(2,0,0);

    StencilWalk query(mesh);
    bool found = query.donorSearch(query_point, start_cell, donor_cell, trace_path);
    if(found)
        std::cout<<"donor cell id: "<<donor_cell<<std::endl;
    tec_out.blocks_ptr.push_back(&mesh);
    tec_out.writeBlocks("mesh.plt");
    tec_out.writeLineSegment("lines.plt", trace_path);
    return 0;
}






