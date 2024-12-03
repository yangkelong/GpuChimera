#include "DonorSearch.h"





template<typename T>
bool StencilWalk<T>::donorSearch(Point query_point, uint32 start_cell, uint32 &donor_cell, std::vector<uint32> &trace_path){
    uint32 cur_cell = start_cell;
    
    //
    Line line;
    uint32 step = 0;
    facet_array = 
    while(step < max_step_){
        trace_path.push_back(cur_cell);

        // 遍历 facet_array, 
        bool any_hit = false;
        for(facet: facet_array){
            if(line.intersect(facet)){
                cur_cell = adj_cell;
                any_hit = true;
                // 取出 adj_cell 单元上的面 facet_array, 并排除刚才击中的面单元 hit_facet
                // 同时 击中的面单元如果为四边形分裂 生成的三角形面元 应该把它的 twin_facet 也排除
                facet_array = 
                break;
            }
        }
        // 如果线段与 facet_array 中的面均不相交 ===> donor_cell = cur_cell
        if(!any_hit){
            donor_cell = cur_cell;
            return true;
        }
        step += 1;
    }
    return false;
}






