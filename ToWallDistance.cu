// ======================================================================== //
// Copyright 2018-2023 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //
#pragma once
#include "cukd/builder.h"
#include "cukd/fcp.h"  // fcp = "find closest point" query
#include <queue>
#include <iomanip>
#include "ToWallDistance.h"
#include <numeric>

using mydata3 = float3;
//using mydata3 = double3;
using mydata = typename cukd::scalar_type_of<mydata3>::type;

struct MyPoint { 
  mydata3  position;
  // 1 byte for split dimension
  uint8_t split_dim; 
};

struct MyPoint_traits : public cukd::default_data_traits<mydata3>{
  using point_t = mydata3;

  static inline __device__ __host__
  mydata3 get_point(const MyPoint &data)
  { return data.position; }

  static inline __device__ __host__
  mydata  get_coord(const MyPoint &data, int dim)
  { return cukd::get_coord(get_point(data),dim); }

  enum { has_explicit_dim = true };

  static inline __device__ void set_dim(MyPoint &p, int dim){p.split_dim = dim; }

  static inline __device__ int  get_dim(const MyPoint &p) { return p.split_dim; }
};

__global__ void d_fcp_mypoint(mydata *d_results, mydata3 *d_queries, int numQueries,
                      const cukd::box_t<mydata3> *d_bounds, MyPoint *d_nodes,
                      int numNodes, mydata cutOffRadius, int *d_records){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= numQueries)
    return;
  mydata3 queryPos = d_queries[tid];  // 查询点 坐标
  cukd::FcpSearchParams params;
  params.cutOffRadius = cutOffRadius;
  int traverse_node_num = 0;
  int closestID = cukd::cct::fcp<MyPoint, MyPoint_traits>(queryPos, *d_bounds, d_nodes, numNodes, params, d_records[tid]);
  // d_records[tid] = traverse_node_num;
  d_results[tid] = (closestID < 0)
                       ? INFINITY
                       : cukd::distance(queryPos, d_nodes[closestID].position);
}


MyPoint* uploadMyPoints(Point *coords, int N) {
    MyPoint* d_points = 0;
    cudaMallocManaged((char**)&d_points, N * sizeof(*d_points));
    if (!d_points)
        throw std::runtime_error("could not allocate points mem...");
    for (int i = 0; i < N; i++) {
        d_points[i].position.x = (mydata) coords[i].x;
        d_points[i].position.y = (mydata) coords[i].y;
        d_points[i].position.z = (mydata) coords[i].z;
    }
    return d_points;
}


template<typename T3>
T3* uploadPoints(Point *coords, int N) {
    T3* d_points = 0;
    cudaMallocManaged((char**)&d_points, N * sizeof(*d_points));
    if (!d_points)
        throw std::runtime_error("could not allocate points mem...");
    for (int i = 0; i < N; i++) {
        d_points[i].x = (mydata) coords[i].x;
        d_points[i].y = (mydata) coords[i].y;
        d_points[i].z = (mydata) coords[i].z;
    }
    return d_points;
}

// 查询host函数, cct
__global__ void d_fcp(mydata *d_results, mydata3 *d_queries, int numQueries,
                      /*! the world bounding box computed by the builder */
                      const cukd::box_t<mydata3> *d_bounds, mydata3 *d_nodes,
                      int numNodes, mydata cutOffRadius, int *d_records){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= numQueries)
    return;
  mydata3 queryPos = d_queries[tid];  // 查询点 坐标
  cukd::FcpSearchParams params;
  params.cutOffRadius = cutOffRadius;
  int traverse_node_num = 0;
  int closestID = cukd::cct::fcp(queryPos, *d_bounds, d_nodes, numNodes, params, d_records[tid]);
  // d_records[tid] = traverse_node_num;
  d_results[tid] = (closestID < 0)
                       ? INFINITY
                       : cukd::distance(queryPos, d_nodes[closestID]);
}

__global__ void d_fcp_stackBased(mydata *d_results, mydata3 *d_queries, int numQueries,
                      /*! the world bounding box computed by the builder */
                      const cukd::box_t<mydata3> *d_bounds, mydata3 *d_nodes,
                      int numNodes, mydata cutOffRadius, int *d_records){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= numQueries)
    return;
  mydata3 queryPos = d_queries[tid];  // 查询点 坐标
  cukd::FcpSearchParams params;
  params.cutOffRadius = cutOffRadius;
  int traverse_node_num = 0;
  int closestID = cukd::stackBased::fcp(queryPos, d_nodes, numNodes, params, d_records[tid]);
  // d_records[tid] = traverse_node_num;
  d_results[tid] = (closestID < 0)
                       ? INFINITY
                       : cukd::distance(queryPos, d_nodes[closestID]);
}


// 后续考虑将构建 查询 分开？
void calWalDist();
extern "C"  void toWallDistance(Point* coords, unsigned int n, Point* query_coords, unsigned int numQueries,
    double* result_coords) {
    std::cout << "toWallDistance ..." <<"k-d tree node_num: "<< n<<", query points: "<< numQueries << std::endl;
    using namespace cukd::common;
    cukd::box_t<mydata3>* d_bounds;
    mydata3* d_points = uploadPoints<mydata3>(coords, n);
    cudaMallocManaged((void**)&d_bounds, sizeof(cukd::box_t<mydata3>));
    std::cout << "allocated memory for the world space bounding box ..." << std::endl;
    // ==================================================================
    // build the tree.
    // ==================================================================
    std::cout << "calling builder..." << std::endl;
    double t0 = cukd::common::getCurrentTime();
    cukd::buildTree(d_points, n, d_bounds);
    CUKD_CUDA_SYNC_CHECK();
    double t1 = cukd::common::getCurrentTime();
    std::cout << "done building tree, took "
        << cukd::common::prettyDouble(t1 - t0) << "s" << std::endl;
    // 搜索时的最大半径
    mydata cutOffRadius = std::numeric_limits<mydata>::infinity();
    mydata3* d_queries = uploadPoints<mydata3>(query_coords, numQueries);
    // allocate memory for the results
    mydata* d_results;
    CUKD_CUDA_CALL(MallocManaged((void**)&d_results, numQueries * sizeof(*d_results)));
    // ==================================================================
    // do queryies 
    // ==================================================================
    // 记录每个查询点 遍历 tree 时 访问的节点数目
    int *d_records;
    cudaMallocManaged((char **)&d_records, numQueries * sizeof(int));
    {
        for (int i = 0; i < numQueries; i++){
            d_records[i] = 0;
        }
        t0 = cukd::common::getCurrentTime();
        int bs = 128;
        int nb = cukd::divRoundUp((int)numQueries, bs);
        d_fcp << <nb, bs >> > (d_results, d_queries, numQueries, d_bounds, d_points, n, cutOffRadius, d_records);
        cudaDeviceSynchronize();
        CUKD_CUDA_SYNC_CHECK();
        t1 = cukd::common::getCurrentTime();
        std::cout << "done "
            << " iterations of " << numQueries
            << " fcp queries, took " << cukd::common::prettyDouble(t1 - t0)
            << "s" << std::endl;
        std::cout << "that is " << cukd::common::prettyDouble(numQueries / (t1 - t0))
            << " queries/s" << std::endl;
        double avg_per_query = std::accumulate(d_records, d_records+numQueries, 0.)/numQueries;
        std::cout << "average traverse_node_num per query: " << avg_per_query << std::endl; 
        for (int i = 0; i < numQueries; i++) {
            result_coords[i] = (double) d_results[i];
        }
    }
    {
        for (int i = 0; i < numQueries; i++){
            d_records[i] = 0;
        }
        t0 = cukd::common::getCurrentTime();
        int bs = 128;
        int nb = cukd::divRoundUp((int)numQueries, bs);
        d_fcp_stackBased << <nb, bs >> > (d_results, d_queries, numQueries, d_bounds, d_points, n, cutOffRadius, d_records);
        cudaDeviceSynchronize();
        CUKD_CUDA_SYNC_CHECK();
        t1 = cukd::common::getCurrentTime();
        std::cout << "done "
            << " iterations of " << numQueries
            << " fcp queries, took " << cukd::common::prettyDouble(t1 - t0)
            << "s" << std::endl;
        std::cout << "that is " << cukd::common::prettyDouble(numQueries / (t1 - t0))
            << " queries/s" << std::endl;
        double avg_per_query = std::accumulate(d_records, d_records+numQueries, 0.)/numQueries;
        std::cout << "average traverse_node_num per query: " << avg_per_query << std::endl; 
        for (int i = 0; i < numQueries; i++) {
            result_coords[i] = (double) d_results[i];
        }
    }

  { 
    //
    cukd::box_t<mydata3> *d_bounds;
    cudaMallocManaged((void **)&d_bounds, sizeof(cukd::box_t<mydata3>));
    MyPoint *d_points = uploadMyPoints(coords, n);
    cukd::buildTree<MyPoint, MyPoint_traits>(d_points, n, d_bounds);
    for (int i = 0; i < numQueries; i++){
      d_records[i] = 0;
    }
    t0 = cukd::common::getCurrentTime();
    int bs = 128;
    int nb = cukd::divRoundUp((int)numQueries, bs);
    d_fcp_mypoint << <nb, bs >> > (d_results, d_queries, numQueries, d_bounds, d_points, n, cutOffRadius, d_records);
    cudaDeviceSynchronize();
    CUKD_CUDA_SYNC_CHECK();
    t1 = cukd::common::getCurrentTime();
    std::cout << "done "
        << " iterations of " << numQueries
        << " fcp queries, took " << cukd::common::prettyDouble(t1 - t0)
        << "s" << std::endl;
    std::cout << "that is " << cukd::common::prettyDouble(numQueries / (t1 - t0))
        << " queries/s" << std::endl;
    double avg_per_query = std::accumulate(d_records, d_records+numQueries, 0.)/numQueries;
    std::cout << "average traverse_node_num per query: " << avg_per_query << std::endl; 
    for (int i = 0; i < numQueries; i++) {
        result_coords[i] = (double) d_results[i];
    }
  }

    cudaFree(d_points);
    cudaFree(d_bounds);
    cudaFree(d_queries);
    cudaFree(d_results);
}

