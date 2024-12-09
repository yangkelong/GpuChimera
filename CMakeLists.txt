# ======================================================================== #
# FFFFFFFFFFFF
# FFFFFFFFFFFF
# FFFFFFFFFFFF
# FFFFFFFFFFFF
# ======================================================================== #
# 额外三方库 gdt cgal boost tecio eigen

cmake_minimum_required(VERSION 3.15...3.31)

project(
  GPUChimera
  VERSION 1.0
  LANGUAGES C CXX CUDA)

set(CMAKE_CUDA_ARCHITECTURES 70)


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ------------------------------------------------------------------
# first, include gdt project to do some general configuration stuff
# ------------------------------------------------------------------
set(gdt_dir ${PROJECT_SOURCE_DIR}/external/gdt/)
# ------------------------------------------------------------------
# import gdt submodule
# ------------------------------------------------------------------
include_directories(${gdt_dir})
# EXCLUDE_FROM_ALL 子目录下的子项目不会被 build
add_subdirectory(${gdt_dir} EXCLUDE_FROM_ALL)
# Eigen
include_directories(${PROJECT_SOURCE_DIR}/external/eigen-3.4.0)
# boost("header only library", no need link)
include_directories(${PROJECT_SOURCE_DIR}/external/boost_1_78_0)
# CGAL("header only library", no need link)
# linux 下 cgal 头文件放入 /usr/local/include/CGAL/~  编译器会默认加入该搜索路径
# 或者都使用项目路径下的 ${PROJECT_SOURCE_DIR}/external/cgal/CGAL
set(cgal_dir ${PROJECT_SOURCE_DIR}/external/cgal-5.6.1/)
include_directories(${cgal_dir}/include)


include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/IO)
# 添加可执行文件
add_executable(my_target 
  GridBlock.h
  DonorSearch.h
  main.cpp
  GridBlock.cpp
  IO/TecplotIO.cpp
  DonorSearch.cpp
  ToWallDistance.h 
  ToWallDistance.cu
  )

# 设置 TECIO 库的安装目录
# 在 Windows 系统中，路径分隔符应使用正斜杠 / 或双反斜杠 \\，而不是单反斜杠 \
# set(TECIO_DIR "C:/Program Files/Tecplot/Tecplot 360 EX 2019 R1/")
set(TECIO_DIR ${PROJECT_SOURCE_DIR}/external/tecio/)
# 查找 TECIO 库
# 注意：这里假设库文件在 TECIO_DIR/lib 目录下，并且库文件名为 tecio.lib（Windows）或 libtecio.a（Linux）
find_library(tecio_LIBRARY NAMES tecio PATHS ${TECIO_DIR}/lib NO_DEFAULT_PATH)
# 为目标添加包含目录
target_include_directories(my_target PUBLIC ${TECIO_DIR}/include)
# 链接 TECIO 库
target_link_libraries(my_target PUBLIC 
	${tecio_LIBRARY} 
	gdt
	)
# windows下 cgal 需要一些额外库: gmp
if(WIN32)
    # gmp
    # auxiliary\gmp\lib
    set(gmp_dir ${cgal_dir}/auxiliary/gmp/)
    include_directories(${gmp_dir}/include)
    target_link_libraries(my_target PUBLIC 
      ${gmp_dir}/lib/libgmp-10.lib
      ${gmp_dir}/lib/libmpfr-4.lib
      )
endif()

# gpu_wall_tree
add_subdirectory(${PROJECT_SOURCE_DIR}/GpuUtils/cuda_kdtree)
target_link_libraries(my_target PUBLIC cudaKDTree)

add_executable(wall_dist ToWallDistance.h ToWallDistance.cu ToWallDistance_main.cpp)
target_link_libraries(wall_dist cudaKDTree)
