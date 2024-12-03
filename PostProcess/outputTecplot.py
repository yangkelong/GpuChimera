import time
from math import pi, cos, sin, sqrt, fabs, ceil
from os import getcwd
import numpy as np
from numba import cuda
import pickle
from multiprocessing import Pool
from DubheGpuSolver.setting.profile import numbaCache
import re

# 函数采用驼峰写法 , 变量采用下划线写法 , 全局变量大写
# 约定 连通性矩阵中存储的编号都是从1开始，有关存储node face cell 连通矩阵中的值均为编号
# 计算数值精度
precision = np.float64
intType = np.int32

# 输出文件目录
OUTPUT_DIR = "output"

# 网格单元维度 DIM = 3 三角形  DIM = 4 长方形
DIM = 3
DELIMITER = "\\"
GAMMA = 1.4
CFL = 2
K2 = 0.75
K4 = 0.02
MA = 0.0
ALPHA = 0.0
R = 287.053
C = 340
ITERMAX = 1000
T = 300
ERROR = 1e-6


# ～计算升力系数和阻力系数
@cuda.jit
def calculateClCd(new_w, faces_dsnxnynz, faces_wall, faces_cell, fxfyfz, n):
    x = cuda.grid(1)
    #
    if x < n:
        no_face = faces_wall[x]
        n = no_face - 1
        ds, nx, ny, nz = faces_dsnxnynz[n, 0], faces_dsnxnynz[n, 1], faces_dsnxnynz[n, 2], faces_dsnxnynz[n, 3]
        fxfyfz[x, 0] = 0.0
        fxfyfz[x, 1] = 0.0
        fxfyfz[x, 2] = 0.0
        no_leftcell = faces_cell[n, 0]
        w = new_w[no_leftcell - 1]
        rho = w[0]
        vx = w[1] / rho
        vy = w[2] / rho
        vz = w[3] / rho
        p = (w[4] - 0.5 * rho * (pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0))) * (GAMMA - 1)
        f = p * ds
        fx = f * nx
        fy = f * ny
        fz = f * nz
        fxfyfz[x, 0] = fx
        fxfyfz[x, 1] = fy
        fxfyfz[x, 2] = fz


#
def output_contour(w_new, XYZ, vertices_amount, cells_node, cell_num, fname_contour='contour.plt'):
    with open(fname_contour, "w") as fname_Con:
        fname_Con.writelines("TITLE =\"3D FLOW\"" + "\n")
        fname_Con.writelines(
            "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Vx\",\"Vy\",\"Vz\",\"V\",\"M\",\"Pressure\",\"Temperature\"" + "\n")
        fname_Con.writelines(
            "ZONE NODES=%d, ELEMENTS=%d, ZONETYPE=FETETRAHEDRON, DATAPACKING=POINT" % (vertices_amount, cell_num) + "\n")
        fname_Con.writelines("DATAPACKING=BLOCK,VARLOCATION=([1-3]=NODAL,[4-11]=CELLCENTERED)" + "\n")
        tmpMatrix = np.zeros((cell_num, 7), dtype=precision)
        for i in range(vertices_amount):
            fname_Con.writelines("%10.5f" % XYZ[i][0] + "\n")
        for i in range(vertices_amount):
            fname_Con.writelines("%10.5f" % XYZ[i][1] + "\n")
        for i in range(vertices_amount):
            fname_Con.writelines("%10.5f" % XYZ[i][2] + "\n")
        for i in range(cell_num):
            rho_f = w_new[i, 0]
            vx_f = w_new[i, 1] / rho_f
            vy_f = w_new[i, 2] / rho_f
            vz_f = w_new[i, 3] / rho_f
            v2 = pow(vx_f, 2) + pow(vy_f, 2) + pow(vz_f, 2)
            v = sqrt(v2)
            p_f = (w_new[i, 4] - 0.5 * rho_f * v2) * (GAMMA - 1)
            Ma_f = v / sqrt(GAMMA * p_f / rho_f)
            tmpMatrix[i, 0] = rho_f
            tmpMatrix[i, 1] = vx_f
            tmpMatrix[i, 2] = vy_f
            tmpMatrix[i, 3] = vz_f
            tmpMatrix[i, 4] = p_f
            tmpMatrix[i, 5] = Ma_f
            tmpMatrix[i, 6] = v
            fname_Con.writelines("%10.5f" % tmpMatrix[i, 0] + "\n")
        for i in range(cell_num):
            fname_Con.writelines("%10.5f" % tmpMatrix[i, 1] + "\n")
        for i in range(cell_num):
            fname_Con.writelines("%10.5f" % tmpMatrix[i, 2] + "\n")
        for i in range(cell_num):
            fname_Con.writelines("%10.5f" % tmpMatrix[i, 3] + "\n")
        for i in range(cell_num):
            fname_Con.writelines("%15.5f" % tmpMatrix[i, 6] + "\n")
        for i in range(cell_num):
            fname_Con.writelines("%15.5f" % tmpMatrix[i, 5] + "\n")
        for i in range(cell_num):
            fname_Con.writelines("%15.5f" % tmpMatrix[i, 4] + "\n")
        for i in range(cell_num):
            T_f = tmpMatrix[i, 4] / (tmpMatrix[i, 0] * R)
            fname_Con.writelines("%15.5f" % T_f + "\n")
        for i in range(cell_num):
            point0, point1, point2, point3 = cells_node[i][0], cells_node[i][1], cells_node[i][2], cells_node[i][3]
            fname_Con.writelines("%6d" % point0 + "\n")
            fname_Con.writelines("%6d" % point1 + "\n")
            fname_Con.writelines("%6d" % point2 + "\n")
            fname_Con.writelines("%6d" % point3 + "\n")


def preOutput(w_new):
    s = time.perf_counter()
    cells_amount = w_new.shape[0]
    threads_per_alock = 512
    blocks_per_grid = ceil(cells_amount / threads_per_alock)
    mtx = np.zeros((cells_amount, 8), dtype=precision)
    w_d = cuda.to_device(w_new)
    mtx_d = cuda.to_device(mtx)
    reget[blocks_per_grid, threads_per_alock](w_d, mtx_d, cells_amount)
    cuda.synchronize()
    tmpMatrix = mtx_d.copy_to_host()
    print("GPU process done! spend time: %f" % (time.perf_counter() - s))
    return tmpMatrix


gamma_Minus_One = GAMMA-1


# 根据守恒量计算原参量等信息
@cuda.jit
def reget(n, cells_W, tmpMatrix, crt):
    i = cuda.grid(1)
    if i < n:
        cell_no = i+1
        convar = cells_W[cell_no]
        rho = convar[0]
        u, v, w = convar[1] / convar[0], convar[2] / convar[0], convar[3] / convar[0]
        W = convar
        p = gamma_Minus_One * (W[4] - 0.5 * (W[1] * W[1] + W[2] * W[2] + W[3] * W[3]) / W[0])
        rho_f = crt[0] * rho
        vx_f = crt[1] * u
        vy_f = crt[1] * v
        vz_f = crt[1] * w
        v2 = pow(vx_f, 2) + pow(vy_f, 2) + pow(vz_f, 2)
        p_f = crt[2] * p
        v = sqrt(v2)
        Ma_f = v / sqrt(GAMMA * p_f / rho_f)
        T_f = p_f / (rho_f * R)
        tmpMatrix[cell_no, 0] = rho_f
        tmpMatrix[cell_no, 1] = vx_f
        tmpMatrix[cell_no, 2] = vy_f
        tmpMatrix[cell_no, 3] = vz_f
        tmpMatrix[cell_no, 4] = v
        tmpMatrix[cell_no, 5] = Ma_f
        tmpMatrix[cell_no, 6] = p_f
        tmpMatrix[cell_no, 7] = T_f


def output_contour(w_new, XYZ, vertices_amount, cells_node, cells_amount, crt, fname_contour='contour.plt'):
    # with open('grid/w.data', 'wb') as f:
    #     f.write(pickle.dumps((w_new, XYZ, vertices_amount, cells_node, cells_amount)))
    # return
    # with open('grid/grid.data', 'rb') as f:
    #     w_new, XYZ, vertices_amount, cells_node, cells_amount = pickle.loads(f.read())
    s = time.perf_counter()
    threads_per_alock = 512
    blocks_per_grid = ceil(cells_amount / threads_per_alock)
    w_d = w_new
    mtx_d = cuda.to_device(np.zeros((cells_amount+1, 8), dtype=precision))
    # 守恒量转化为需要的变量
    crt_d = cuda.to_device(crt)
    reget[blocks_per_grid, threads_per_alock](cells_amount, w_d, mtx_d, crt_d)
    cuda.synchronize()
    tmpMatrix = mtx_d.copy_to_host()
    with open(fname_contour, "w") as fname_Con:
        fname_Con.writelines("TITLE =\"3D FLOW\"" + "\n")
        fname_Con.writelines(
            "VARIABLES=\"X\",\"Y\",\"Z\",\"Density\",\"Vx\",\"Vy\",\"Vz\",\"V\",\"M\",\"Pressure\",\"Temperature\"" + "\n")
        fname_Con.writelines("ZONE NODES=%d, ELEMENTS=%d, ZONETYPE=FETETRAHEDRON" % (vertices_amount, cells_amount) + "\n")
        fname_Con.writelines("DATAPACKING=BLOCK,VARLOCATION=([1-3]=NODAL,[4-11]=CELLCENTERED)" + "\n")
        x_block = ''
        y_block = ''
        z_block = ''
        rho_block = ''
        vx_block = ''
        vy_block = ''
        vz_block = ''
        v_block = ''
        m_block = ''
        p_block = ''
        t_block = ''
        points_block = ''
        for i in range(1, vertices_amount+1):
            x_block += "%10.5f" % XYZ[i][0] + "\n"
            y_block += "%10.5f" % XYZ[i][1] + "\n"
            z_block += "%10.5f" % XYZ[i][2] + "\n"
            # x_block += f"{XYZ[i][0]:15.5f}" + "\n"
            # y_block += f"{XYZ[i][1]:15.5f}" + "\n"
            # z_block += f"{XYZ[i][2]:15.5f}" + "\n"
        for i in range(1, cells_amount+1):
            rho_block += "%15.5f" % tmpMatrix[i, 0] + "\n"
            vx_block += "%15.5f" % tmpMatrix[i, 1] + "\n"
            vy_block += "%15.5f" % tmpMatrix[i, 2] + "\n"
            vz_block += "%15.5f" % tmpMatrix[i, 3] + "\n"
            v_block += "%15.5f" % tmpMatrix[i, 4] + "\n"
            m_block += "%15.5f" % tmpMatrix[i, 5] + "\n"
            p_block += "%15.5f" % tmpMatrix[i, 6] + "\n"
            t_block += "%15.5f" % tmpMatrix[i, 7] + "\n"
            # rho_block += f"{tmpMatrix[i, 0]:15.5f}" + "\n"
            # vx_block += f"{tmpMatrix[i, 1]:15.5f}" + "\n"
            # vy_block += f"{tmpMatrix[i, 2]:15.5f}" + "\n"
            # vz_block += f"{tmpMatrix[i, 3]:15.5f}" + "\n"
            # v_block += f"{tmpMatrix[i, 4]:15.5f}" + "\n"
            # m_block += f"{tmpMatrix[i, 5]:15.5f}" + "\n"
            # p_block += f"{tmpMatrix[i, 6]:15.5f}" + "\n"
            # t_block += f"{tmpMatrix[i, 7]:15.5f}" + "\n"
            cells_node[i] = list(cells_node[i])
            point0, point1, point2, point3 = cells_node[i][0], cells_node[i][1], cells_node[i][2], cells_node[i][3]
            points_block += f"{point0:6d}, {point1:6d}, {point2:6d}, {point3:6d}" + "\n"
        blocks = [x_block, y_block, z_block, rho_block, vx_block, vy_block, vz_block, v_block, m_block, p_block,
                  t_block, points_block]
        for b in blocks:
            fname_Con.write(b)
        print("write done! spend time: %f" % (time.perf_counter() - s))



# 测试函数，将右半侧 x>0 区域守恒量赋值，测试tecplot输出是否有问题
@cuda.jit(cache=numbaCache)
def setConvar(n, cells_W, cells_centerCoord):
    x = cuda.grid(1)
    if x < n:
        cell_no = x+1
        if cells_centerCoord[cell_no, 0] > 0.:
            cells_W[cell_no] = (1.0, 0.0, 0.0, 0.0, 1.7857143)


if __name__ == "__main__":
    s = time.perf_counter()
    print("spend time:%f" % (time.perf_counter() - s))
