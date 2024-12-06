#include "DonorSearch.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;

bool StencilWalk::donorSearch(const Point &query_point, uint32 start_cell, uint32 &donor_cell, std::vector<Point> &trace_path){
    uint32 cur_cell = start_cell;
    // 定义线段
    Point p_start = mesh.cells_center[start_cell];
    Point_3 start_point(p_start.x, p_start.y, p_start.z);
    Point_3 end_point(query_point.x, query_point.y, query_point.z);
    Segment_3 seg(start_point, end_point);
    uint32 step = 0;
    std::set<uint32> facet_set;
    mesh.getCellTriFacet(cur_cell, facet_set);
    uint32 intersection_count = 0;  // 线段与三角形相交判断计数
    while(step < max_step){
        trace_path.push_back(mesh.cells_center[cur_cell]);
        gdt::printVec(mesh.cells_center[cur_cell]);
        // 遍历 facet_set, 
        bool any_hit = false;
        for(uint32 facet_id: facet_set){
            // 能不能直接利用 vec3f 定义 Point_3？？？？？
            const auto &tri_facet = mesh.tri_facet_array[facet_id];
            Point p_0 = mesh.vertex_coord[tri_facet.index.x];
            Point p_1 = mesh.vertex_coord[tri_facet.index.y];
            Point p_2 = mesh.vertex_coord[tri_facet.index.z];
            // 定义点
            Point_3 p1(p_0.x, p_0.y, p_0.z);
            Point_3 p2(p_1.x, p_1.y, p_1.z);
            Point_3 p3(p_2.x, p_2.y, p_2.z);
            // 定义三角形
            Triangle_3 tri(p1, p2, p3);
            // 求交
            CGAL::cpp11::result_of<K::Intersect_3(Segment_3, Triangle_3)>::type result = CGAL::intersection(seg, tri);
            intersection_count += 1;
            if(result){
                /*  交点如果刚好位于边上 有问题...
                    如果交点位于边上 则算出多个 adj_cells 然后在这些 adj_cells 中选一个最合适的
                    例如: adj_cell_center 离 query_point 最近 or 离 seg 最近
                */ 
                std::cout << "Intersection point: " << *boost::get<Point_3>(&*result) << std::endl;
                const auto& side_cells = mesh.face_side_cell[tri_facet.parent_facet_id];
                const uint32 adj_cell = side_cells.x==cur_cell ? side_cells.y : side_cells.x;
                cur_cell = adj_cell;
                any_hit = true;
                // 取出 adj_cell 单元上的面 facet_set, 并排除刚才击中的面单元 hit_facet
                // 同时 击中的面单元如果为四边形分裂 生成的三角形面元 应该把它的 twin_facet 也排除
                mesh.getCellTriFacet(cur_cell, facet_set);
                facet_set.erase(facet_id);
                if(tri_facet.have_twin) facet_set.erase(tri_facet.twin_facet_id);
                break;
            }
        }
        // 如果线段与 facet_set 中的面均不相交 ===> donor_cell = cur_cell
        if(!any_hit){
            donor_cell = cur_cell;
            return true;
        }
        step += 1;
    }
    return false;
}






