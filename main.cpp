#include "GridBlock.h"
#include "IO/TecplotIO.h"

// #include "TecplotWriter.h"


int main(int argc, char *argv[]) {
    Block grid(argv[1]);
    TecplotIO tec_out;
    tec_out.blocks_ptr.push_back(&grid);
    tec_out.writeBlocks("mesh.plt");
    return 0;
}






