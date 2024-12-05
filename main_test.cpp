#include "GridBlock.h"
// #include "TecplotWriter.h"
#include "TECIO.h"
#include "MASTER.h" /* for defintion of NULL */
#include "gdt/math/vec.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;

#include <Eigen/Dense>

int main(int argc, char *argv[]) {
    gdt::vec3f  point_1{0,1,2};
    gdt::vec3f  point_2{3,4,5};
    gdt::vec3f  point_3 = gdt::cross(point_1, point_2);
    gdt::printVec(point_3);


    // 定义点
    Point_3 p1(0, 0, 0);
    Point_3 p2(1, 0, 0);
    Point_3 p3(0, 1, 0);
    Point_3 start_point(0.15, 0.25, -1);
    Point_3 end_point(0.25, 0.15, 1);
    // 定义线段和三角形
    Segment_3 seg(start_point, end_point);
    Triangle_3 tri(p1, p2, p3);

    // 求交
    CGAL::cpp11::result_of<K::Intersect_3(Segment_3, Triangle_3)>::type
        result = CGAL::intersection(seg, tri);

    // 检查交点并输出
    if (result) {
        if (const Point_3* p = boost::get<Point_3>(&*result)) {
            std::cout << "Intersection point: " << *p << std::endl;
        }
        else if (const Segment_3* s = boost::get<Segment_3>(&*result)) {
            std::cout << "Intersection segment: " << *s << std::endl;
        }
    }
    else {
        std::cout << "No intersection" << std::endl;
    }
    using Eigen::Matrix;
    using Eigen::Vector;
    // 创建一个 2x2 的矩阵
    Matrix<int, 2, 2> mat;
    mat << 1, 2, 3, 4;
    std::cout << "Matrix:\n" << mat << std::endl;

    // 创建一个 2D 向量
    Vector<float, 2> vec(99.9, 88.8);
    std::cout << "Vector:\n" << vec << std::endl;

    // 动态大小的矩阵，使用列表初始化
    Eigen::MatrixXi  dymat(2, 3);
    dymat << 1, 2, 3, 4, 5, 6;
    std::cout << "Dynamic matrix:\n" << dymat << std::endl;

    // 矩阵乘法
    Eigen::VectorXi res = mat * Eigen::VectorXi::LinSpaced(2, 1, 2);
    std::cout << "Matrix-vector multiplication result:\n" << res << std::endl;

    return 0;
    Block grid(argv[1]);
    /* Call TECINI142 */
    INTEGER4 FileType = 0;   /* 0 for full file */
    INTEGER4 FileFormat = 0; // 0 == PLT, 1 == SZPLT; Only PLT is currently
    INTEGER4 Debug = 0;
    INTEGER4 VIsDouble = 1;
    INTEGER4 res_code = 0;                      /* use to check return codes */
    // 初始化
    res_code = TECINI142((char *)"Pyramid",     /* Data Set Title */
                  (char *)"X Y Z",       /* Variable List */
                  (char *)"pyramid.plt", /* File Name */
                  (char *)".",           /* Scratch Directory */
                  &FileFormat,
                  &(FileType),
                  &(Debug),
                  &(VIsDouble));
    /* Call TECZNE142 */
    INTEGER4 ZoneType = 7;   /* 7 for FEPolyhedron */
    INTEGER4 NumNodes = 5;   /* number of unique nodes */
    INTEGER4 NumElems = 1;   /* number of elements */
    INTEGER4 NumFaces = 5;   /* number of unique faces */
    INTEGER4 ICellMax = 0;   /* Not Used, set to zero */
    INTEGER4 JCellMax = 0;   /* Not Used, set to zero */
    INTEGER4 KCellMax = 0;   /* Not Used, set to zero */
    double SolTime = 12.65;  /* solution time */
    INTEGER4 StrandID = 0;   /* static zone */
    INTEGER4 ParentZone = 0; /* no parent zone */
    INTEGER4 IsBlock = 1;    /* block format */
    INTEGER4 NFConns = 0;    /* not used for FEPolyhedron
                              * zones
                              */
    INTEGER4 FNMode = 0;     /* not used for FEPolyhedron
                              * zones
                              */
    INTEGER4 *PassiveVarArray = NULL;
    INTEGER4 *ValueLocArray = NULL;
    INTEGER4 *VarShareArray = NULL;
    INTEGER4 ShrConn = 0;
    /* The number of face nodes in the zone. This example creates
     * a zone with a single pyramidal cell. This cell has four
     * triangular faces and one rectangular face, yielding a total
     * of 16 face nodes.
     */
    INTEGER4 NumFaceNodes = 16;
    INTEGER4 NumBConns = 0; /* No Boundary Connections */
    INTEGER4 NumBItems = 0; /* No Boundary Items */
    // 创建zone
    res_code = TECZNE142((char *)"Polyhedral Zone (Octahedron)", &ZoneType, &NumNodes, &NumElems, &NumFaces,
                  &ICellMax, &JCellMax, &KCellMax, &SolTime, &StrandID, &ParentZone, &IsBlock,
                  &NFConns, &FNMode, &NumFaceNodes, &NumBConns, &NumBItems, PassiveVarArray,
                  ValueLocArray, VarShareArray, &ShrConn);
    /* Initialize arrays of nodal data */
    double *X = new double[NumNodes];
    double *Y = new double[NumNodes];
    double *Z = new double[NumNodes];
    X[0] = 0;
    Y[0] = 0;
    Z[0] = 0;
    X[1] = 1;
    Y[1] = 1;
    Z[1] = 2;
    X[2] = 2;
    Y[2] = 0;
    Z[2] = 0;
    X[3] = 2;
    Y[3] = 2;
    Z[3] = 0;
    X[4] = 0;
    Y[4] = 2;
    Z[4] = 0;
    /* Write the data (using TECDAT142) */
    INTEGER4 DIsDouble = 1; /* One for double precision */
    res_code = TECDAT142(&NumNodes, X, &DIsDouble);
    res_code = TECDAT142(&NumNodes, Y, &DIsDouble);
    res_code = TECDAT142(&NumNodes, Z, &DIsDouble);
    delete X;
    delete Y;
    delete Z;
    /* Define the Face Nodes.
     *
     * The FaceNodes array is used to indicate which nodes define
     * which face. As mentioned earlier, the number of the nodes is
     * implicitly defined by the order in which the nodal data is
     * provided. The first value of each nodal variable describes
     * node 1, the second value describes node 2, and so on.
     *
     * The face numbering is implicitly defined. Because there are
     * two nodes in each face, the first two nodes provided define
     * face 1, the next two define face 2 and so on. If there was
     * a variable number of nodes used to define the faces, the
     * array would be more complicated.
     */
    INTEGER4 *FaceNodeCounts = new INTEGER4[NumFaces];
    /* The first four faces are triangular, i.e. have three nodes.
     * The fifth face is rectangular, i.e. has four nodes. */
    FaceNodeCounts[0] = 3;
    FaceNodeCounts[1] = 3;
    FaceNodeCounts[2] = 3;
    FaceNodeCounts[3] = 3;
    FaceNodeCounts[4] = 4;
    INTEGER4 *FaceNodes = new INTEGER4[NumFaceNodes];
    /* Face Nodes for Face 1 */
    FaceNodes[0] = 1;
    FaceNodes[1] = 2;
    FaceNodes[2] = 3;
    /* Face Nodes for Face 2 */
    FaceNodes[3] = 3;
    FaceNodes[4] = 2;
    FaceNodes[5] = 4;
    /* Face Nodes for Face 3 */
    FaceNodes[6] = 5;
    FaceNodes[7] = 2;
    FaceNodes[8] = 4;
    /* Face Nodes for Face 4 */
    FaceNodes[9] = 1;
    FaceNodes[10] = 2;
    FaceNodes[11] = 5;
    /* Face Nodes for Face 5 */
    FaceNodes[12] = 1;
    FaceNodes[13] = 5;
    FaceNodes[14] = 4;
    FaceNodes[15] = 3;
    /* Define the right and left elements of each face.
     *
     * The last step for writing out the polyhedral data is to
     * define the right and left neighboring elements for each
     * face. The neighboring elements can be determined using the
     * right-hand rule. For each face, place your right-hand along
     * the face which your fingers pointing the direction of
     * incrementing node numbers (i.e. from node 1 to node 2).
     * Your right thumb will point towards the right element; the
     * element on the other side of your hand is the left element.
     *
     * The number zero is used to indicate that there isn’t an
     * element on that side of the face.
     *
     * Because of the way we numbered the nodes and faces, the
     * right element for every face is the element itself
     * (element 1) and the left element is “no-neighboring element"
     * (element 0).
     */
    INTEGER4 *FaceLeftElems = new INTEGER4[NumFaces];
    FaceLeftElems[0] = 1;
    FaceLeftElems[1] = 1;
    FaceLeftElems[2] = 0;
    FaceLeftElems[3] = 0;
    FaceLeftElems[4] = 0;
    INTEGER4 *FaceRightElems = new INTEGER4[NumFaces];
    FaceRightElems[0] = 0;
    FaceRightElems[1] = 0;
    FaceRightElems[2] = 1;
    FaceRightElems[3] = 1;
    FaceRightElems[4] = 1;
    /* Write the face map (created above) using TECPOLYFACE142. */
    res_code = TECPOLYFACE142(&NumFaces,
                       FaceNodeCounts,  /* The face node counts array */
                       FaceNodes,       /* The face nodes array */
                       FaceLeftElems,   /* The left elements array */
                       FaceRightElems); /* The right elements array */
    delete FaceNodeCounts;
    delete FaceNodes;
    delete FaceLeftElems;
    delete FaceRightElems;
    res_code = TECEND142();
    return 0;
}






