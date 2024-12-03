# 残差, 升力系数, 阻力系数等参数
import numpy as np
from DubheGpuSolver.Communication.MyMPI import MyMPI
from DubheGpuSolver.setting.profile import OutputDir, nVar, nVar_Turb, THREADSPERBLOCK, nDIM, Ref_area, \
    Ref_DynamicPressure, Crt_LENGTH, Crt_PRESSURE, PRECISION, FreestreamDirection, AOA, MPI_Root, MPI_isRoot, \
    Ref_MACH, Re, Ref_area_NonDim, Cgamma, MPI_Rank
import pickle
from DubheGpuSolver.init import cuda, ceil, cupy, sqrt, radians, cos, sin
from DubheGpuSolver.Solver.Solver import getPressure_device_fn, getViscosity_device_fn, getTemperature_device_fn


def writeSolution(cells_W, file="fluid.data"):
    file = OutputDir + file
    with open(file, "wb") as f:
        f.write(pickle.dumps(cells_W))


def loadSolutiion(file="fluid.data"):
    file = OutputDir + file
    with open(file, "rb") as f:
        return pickle.loads(f.read())


# 压强作用力系数 * 动压 * 参考面积 = 压强作用力
@cuda.jit
def calPressureForceCoff(n, faces_no, cells_W, faces_sideCellsNo, faces_normal, faces_area, cells_F):
    x = cuda.grid(1)
    #
    if x < n:
        face_no = faces_no[x]
        UnitNormal = faces_normal[face_no]
        Area = faces_area[face_no]
        c0, c1 = faces_sideCellsNo[face_no]
        W_i = cells_W[c0]
        p = getPressure_device_fn(W_i)
        # 该处是否有必要 - 1/Cgamma
        p = p - 1/Cgamma
        # 计算压强作用在壁面上的力系数
        dFp = p * Area / (0.5*Ref_MACH*Ref_MACH * Ref_area_NonDim)
        for i in range(nDIM):
            cells_F[x, i] = dFp * UnitNormal[i]
        pass


@cuda.jit
def calViscousForceCoff(n, faces_no, cells_W, cells_VelT_Gdt, faces_sideCellsNo, faces_normal, faces_area, cells_F):
    x = cuda.grid(1)
    #
    if x < n:
        face_no = faces_no[x]
        UnitNormal = faces_normal[face_no]
        Area = faces_area[face_no]
        c0, c1 = faces_sideCellsNo[face_no]
        W_i = cells_W[c0]
        # 粘性力
        # /* mju_Turb on wall face is 0.0 */
        T_i = getTemperature_device_fn(W_i)
        mju_Laminar = getViscosity_device_fn(T_i)
        viscosity = Ref_MACH / Re * mju_Laminar
        # //--- divergence of velocity ------//
        VelT_Gdt = cells_VelT_Gdt[c0]
        div_vel = 0.
        for i in range(3):
            div_vel += VelT_Gdt[i, i]
        # 应力张量 tau //-------- stress tensor ------------//
        tau = cuda.local.array((nDIM, nDIM), PRECISION)
        for i in range(nDIM):
            tau[i, i] = 2.0 * viscosity * (VelT_Gdt[i, i] - 1.0 / 3.0 * div_vel)
            for j in range(i + 1, nDIM):
                tau[i, j] = tau[j, i] = viscosity * (VelT_Gdt[i, j] + VelT_Gdt[j, i])
        # 壁面上应力 = 应力张量 点乘 法向(指向流场内侧)
        Stress = cuda.local.array((nDIM,), PRECISION)
        Fv = viscousForce = cuda.local.array((nDIM,), PRECISION)
        # 粘性作用在壁面上的力
        for idim in range(nDIM):
            Stress[idim] = 0.
            for jdim in range(nDIM):
                Stress[idim] -= tau[idim][jdim] * UnitNormal[jdim]
            Fv[idim] = Stress[idim] * Area / (0.5*Ref_MACH*Ref_MACH * Ref_area_NonDim)
        # 将粘性力 计入
        for i in range(nDIM):
            cells_F[x, i] = Fv[i]
        pass


class OutputText:
    def __init__(self, faces_wall, cellSet, solve):
        self.isOutputClCd = solve.isOutputClCd
        if not self.isOutputClCd:
            return
        self.n = faces_wall.shape[0]
        self.faces_wall = cuda.to_device(faces_wall)
        # Note: 力系数
        self.cells_ForceCoff = cuda.device_array((self.n, nDIM), dtype=PRECISION)
        self.cells_ForceCoff_cupy = cupy.asarray(self.cells_ForceCoff)
        self.cellSet = cellSet
        self.solve = solve
        self.clCdSet = list()
        # self.Ref_F = 0.5*Ref_MACH*Ref_MACH * Ref_area_NonDim
        self.Ref_F = Ref_area * Ref_DynamicPressure
        self.AerodynamicCoefficient = np.zeros((2,), dtype=PRECISION)
        self.IsViscous = solve.IsViscous
        pass

    # ～计算升力系数和阻力系数
    def calClCd(self):
        if not self.isOutputClCd:
            return
        cl = 0.
        cd = 0.
        # 有可能该分区没有壁面
        if self.n > 0:
            calPressureForceCoff[ceil(self.n / THREADSPERBLOCK), THREADSPERBLOCK](self.n, self.faces_wall,
                                                                                  self.solve.cells_W,
                                                                                  self.solve.faces_sideCellsNo,
                                                                                  self.solve.faces_normal,
                                                                                  self.solve.faces_area, self.cells_ForceCoff)
            cuda.synchronize()
            # 计算合 力系数
            ForceCoff = self.cells_ForceCoff_cupy.sum(axis=0)
            # print("ForceCoff_p", ForceCoff)
            if self.IsViscous > 0:
                # faces_no, cells_W, cells_VelT_Gdt, faces_sideCellsNo, faces_normal, faces_area, cells_F
                calViscousForceCoff[ceil(self.n / THREADSPERBLOCK), THREADSPERBLOCK](self.n, self.faces_wall,
                                                                                     self.solve.cells_W,
                                                                                     self.solve.cells_VelT_Gdt,
                                                                                     self.solve.faces_sideCellsNo,
                                                                                     self.solve.faces_normal,
                                                                                     self.solve.faces_area, self.cells_ForceCoff)
                cuda.synchronize()
                ForceCoff_v = self.cells_ForceCoff_cupy.sum(axis=0)
                ForceCoff += ForceCoff_v
                # print("ForceCoff_v", ForceCoff_v)
            if nDIM == 3:
                cl = ForceCoff[2] * cos(radians(AOA)) - ForceCoff[0] * sin(radians(AOA))
                cd = ForceCoff[2] * sin(radians(AOA)) + ForceCoff[0] * cos(radians(AOA))
            else:
                pass
        self.AerodynamicCoefficient[0] = cl
        self.AerodynamicCoefficient[1] = cd
        # 通信根进程获得 cl, cd
        AerodynamicCoefficient = MyMPI.reduce(self.AerodynamicCoefficient, MPI_Root, 1)
        if MPI_isRoot:
            cl = AerodynamicCoefficient[0]
            cd = AerodynamicCoefficient[1]
            print("cl: %10.6f, cd: %10.6f" % (cl, cd))
        self.clCdSet.append((cl, cd))

    def writeClCd(self, fileName="clcd.dat"):
        if not self.isOutputClCd:
            return
        fname_res = OutputDir + fileName
        with open(fname_res, "w", encoding="utf-8") as fname_Res:
            fname_Res.write("TITLE = \"AeroCoeff-history\"" + "\n")
            fname_Res.write("VARIABLES=  \"Iter\", \"cl\",\"cd\"" + "\n")
            # fname_Res.writelines("%10s%10s%10s" % ("Iter", "cl", "cd") + "\n")
            VarNum = nVar + nVar_Turb
            for iteration, clcd in enumerate(self.clCdSet):
                fname_Res.writelines(("%10s" + "%10.6f" * 2) % (iteration + 1, *clcd) + "\n")

    def writeResidual(self, ResidualLog, isVicious=0, fileName="res.dat"):
        fname_res = OutputDir + fileName
        with open(fname_res, "w", encoding="utf-8") as fname_Res:
            fname_Res.write("TITLE = \"residual-history\"" + "\n")
            if isVicious <= 1:
                fname_Res.write(
                    "VARIABLES=  \"Iter\", \"Res_rho\", \"Res_rho_u\", \"Res_rho_v\", \"Res_rho_w\", \"Res_rho_e\"" + "\n")
                # fname_Res.writelines("%10s%10s%10s%10s%10s%10s" % ("Iter", "Res_rho", "Res_rho_u", "Res_rho_v",
                #                                                    "Res_rho_w", "Res_rho_e") + "\n")
            elif isVicious == 2:
                fname_Res.write(
                    "VARIABLES=  \"Iter\", \"Res_rho\", \"Res_rho_u\", \"Res_rho_v\", \"Res_rho_w\", \"Res_rho_e\", \"Res_turb_nju\"" + "\n")
                # fname_Res.writelines("%10s%10s%10s%10s%10s%10s%15s" % ("Iter", "Res_rho", "Res_rho_u", "Res_rho_u",
                #                                                        "Res_rho_w", "Res_rho_e", "Res_turb_nju") + "\n")
            VarNum = nVar + nVar_Turb
            for iteration, res_it in enumerate(ResidualLog):
                fname_Res.writelines(("%10s" + "%10.6f" * VarNum) % (iteration + 1, *res_it) + "\n")

    def writeIterTimeHistory(self, iterTimeList, fileName="iterTime.dat"):
        fname_res = OutputDir + fileName
        with open(fname_res, "w", encoding="utf-8") as fname_Res:
            fname_Res.write("TITLE = \"IterTime-history\"" + "\n")
            fname_Res.write("VARIABLES= \"Iter\", \"iterTime\"" + "\n")
            for iteration in range(iterTimeList.shape[0]-1):
                iterTime = iterTimeList[iteration+1]
                fname_Res.writelines(("%10s" + "%10.2f") % (iteration + 1, iterTime) + "\n")

    def writeFuildResidual(self, fluidResidualList, isVicious=0, fileName="fluidRes.dat"):
        fname_res = OutputDir + fileName
        with open(fname_res, "w", encoding="utf-8") as fname_Res:
            fname_Res.write("TITLE = \"residual-history\"" + "\n")
            if isVicious <= 1:
                fname_Res.write(
                    "VARIABLES=  \"Iter\", \"res_rho\", \"res_rho_u\", \"res_rho_v\", \"res_rho_w\", \"res_rho_e\"" + "\n")
            elif isVicious == 2:
                fname_Res.write(
                    "VARIABLES=  \"Iter\", \"res_rho\", \"res_rho_u\", \"res_rho_v\", \"res_rho_w\", \"res_rho_e\", \"res_turb_nju\"" + "\n")
            VarNum = nVar + nVar_Turb
            for iteration in range(fluidResidualList.shape[0] - 1):
                res_it = fluidResidualList[iteration + 1]
                fname_Res.writelines(("%10s" + "%10.6f" * VarNum) % (iteration + 1, *res_it) + "\n")