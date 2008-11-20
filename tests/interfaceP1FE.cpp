#include "geom/builder.h"
#include "levelset/levelset.h"

#include <fstream>

using namespace DROPS;


/// \brief Routine to number unknowns on the vertices surrounding an
/// interface.
///
/// This function allocates memory for the Unknown-indices in system
/// idx on all vertices belonging to tetras between begin and end which
/// are cut by the zero level of lset. A tetra is cut if InterfacePatchCL
/// says so.
///
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is #Unknowns+stride.
void CreateNumbOnInterfaceVertex (const Uint idx, IdxT& counter, Uint stride,
    const MultiGridCL::TriangTetraIteratorCL& begin,
    const MultiGridCL::TriangTetraIteratorCL& end,
    const VecDescCL& ls)
{
    if (stride == 0) return;

    InterfacePatchCL p;
    for (MultiGridCL::TriangTetraIteratorCL it= begin; it != end; ++it) {
        p.Init( *it, ls);
        if (!p.Intersects()) continue;

        for (Uint i= 0; i < NumVertsC; ++i) {
            UnknownHandleCL& u= const_cast<VertexCL*>( it->GetVertex( i))->Unknowns;
            u.Prepare( idx);
            if ( u( idx) == NoIdx) {
                u( idx)= counter;
                counter+= stride;
            }
        }
    }
}

template<class BndDataT>
void CreateNumbOnInterface(Uint level, IdxDescCL& idx, MultiGridCL& mg,
    const VecDescCL& ls)
{
    // set up the index description
    idx.SetTriangLevel( level);
    idx.SetNumUnknowns( 0);

    const Uint idxnum= idx.GetIdx();
    IdxT num_unknowns = idx.NumUnknowns();
    // allocate space for indices; number unknowns in TriangLevel level
    if (idx.NumUnknownsVertex)
        CreateNumbOnInterfaceVertex( idxnum, num_unknowns, idx.NumUnknownsVertex(),
            mg.GetTriangTetraBegin( level), mg.GetTriangTetraEnd( level), ls);
    idx.SetNumUnknowns( num_unknowns);
    if (idx.NumUnknownsEdge() != 0 || idx.NumUnknownsFace() != 0 || idx.NumUnknownsTetra() != 0)
        throw DROPSErrCL( "CreateNumbOnInterface: Only vertex unknowns are implemented\n" );
}


void
TestSingleTetra()
{
    TetraBuilderCL tetra( 0); // unrefined reference tetra
    MultiGridCL mg( tetra);

    LevelsetP2CL lset( mg);

    lset.idx.CreateNumbering( 0, mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Phi.Data= 1.0;

    IdxDescCL ifaceidx( P1_FE);
    ifaceidx.SetTriangLevel( 0);
    IdxT num_unknowns = ifaceidx.NumUnknowns();
    std::cout << "Testing vertex numbering around no interface:" << std::endl;
    CreateNumbOnInterfaceVertex( ifaceidx.GetIdx(), num_unknowns,
        ifaceidx.NumUnknownsVertex(), mg.GetTriangTetraBegin( 0), mg.GetTriangTetraEnd( 0),
        lset.Phi);
    ifaceidx.SetNumUnknowns( num_unknowns);
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);

    std::cout << "Testing vertex numbering interface in 1 tetra:\n" << std::endl;
    lset.Phi.Data[0]= -1.0;
    num_unknowns = ifaceidx.NumUnknowns();
    CreateNumbOnInterfaceVertex( ifaceidx.GetIdx(), num_unknowns,
        ifaceidx.NumUnknownsVertex(), mg.GetTriangTetraBegin( 0), mg.GetTriangTetraEnd( 0),
        lset.Phi);
    ifaceidx.SetNumUnknowns( num_unknowns),
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);
}


double plane(const Point3DCL& p)
{
    return p[0] - 0.41;
}

double plane2(const Point3DCL& p)
{
    return p[0] - 0.4;
}

double sphere_1(const Point3DCL& p)
{
    Point3DCL c( 0.47);
    Point3DCL x= fabs( p - c);

    return x[0] + 1.11*x[1] + 1.111*x[2] - 0.31;
}

double sphere_2(const Point3DCL& p)
{
    Point3DCL c( 0.5);
    Point3DCL x= p - c;

    return x.norm() - 0.31;
}

double sphere_max(const Point3DCL& p)
{
    Point3DCL c( 0.47);
    Point3DCL x= fabs( p - c);

    return std::max( x[0], std::max( x[1], x[2])) - 0.31;
}

void
TestPlaneInCube()
{
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                10, 1, 1);
    MultiGridCL mg( brick);

    LevelsetP2CL lset( mg);

    lset.idx.CreateNumbering( 0, mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Phi.Data= 1.0;

    IdxDescCL ifaceidx( P1_FE);
    ifaceidx.SetTriangLevel( 0);
    IdxT num_unknowns = ifaceidx.NumUnknowns();
    std::cout << "Testing vertex numbering around planar interface:" << std::endl;
    lset.Init( plane);
    CreateNumbOnInterfaceVertex( ifaceidx.GetIdx(), num_unknowns,
        ifaceidx.NumUnknownsVertex(), mg.GetTriangTetraBegin( 0), mg.GetTriangTetraEnd( 0),
        lset.Phi);
    ifaceidx.SetNumUnknowns( num_unknowns);
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);

    num_unknowns = ifaceidx.NumUnknowns();
    std::cout << "Testing vertex numbering around planar interface containing vertices:" << std::endl;
    lset.Init( plane2);
    CreateNumbOnInterfaceVertex( ifaceidx.GetIdx(), num_unknowns,
        ifaceidx.NumUnknownsVertex(), mg.GetTriangTetraBegin( 0), mg.GetTriangTetraEnd( 0),
        lset.Phi);
    ifaceidx.SetNumUnknowns( num_unknowns);
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;
    ifaceidx.DeleteNumbering( mg);
}


void SetupInterfaceMassP1OnTriangle (const LocalP1CL<> p1[4],
    Quad5_2DCL<> q[4], MatrixBuilderCL& M, const IdxT Numb[4],
    const BaryCoordCL triangle[3], double det)
{
    for (int i= 0; i < 4; ++i)
        q[i].assign( p1[i], triangle);

    Quad5_2DCL<> m;
    double tmp;
    for (int i= 0; i < 4; ++i) {
        m= q[i]*q[i];
        M( Numb[i], Numb[i])+= m.quad( det);
        for(int j= 0; j < i; ++j) {
            m= (q[j]*q[i]);
            tmp= m.quad( det);
            M( Numb[i], Numb[j])+= tmp;
            M( Numb[j], Numb[i])+= tmp;
        }
    }
}

void SetupInterfaceMassP1 (const MultiGridCL& MG, MatDescCL* matM, const VecDescCL& ls)
{
    const IdxT num_unks_pr=  matM->RowIdx->NumUnknowns();
    MatrixBuilderCL M( &matM->Data, num_unks_pr,  num_unks_pr);

    const Uint lvl= matM->GetRowLevel();
    IdxT Numb[4];

    double det;
    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> q[4], m;
    
    InterfacePatchCL patch;
    DROPS_FOR_TRIANG_CONST_TETRA( MG, lvl, it) {
        patch.Init( *it, ls);

        GetLocalNumbP1NoBnd( Numb, *it, *matM->RowIdx);
        for (int ch= 0; ch < 8; ++ch) {
            if (!patch.ComputeForChild( ch)) // no patch for this child
                continue;

            det= patch.GetFuncDet();
            SetupInterfaceMassP1OnTriangle( p1, q, M, Numb, &patch.GetBary( 0), det);
            if (patch.IsQuadrilateral()) {
                det*= patch.GetAreaFrac();
                SetupInterfaceMassP1OnTriangle( p1, q, M, Numb, &patch.GetBary( 1), det);
            }
        }
   }
    M.Build();
}


int main ()
{
  try {
    TestSingleTetra();
    TestPlaneInCube();

    std::cerr << "SetupInterfaceMassP1:\n";
    DROPS::BrickBuilderCL brick(DROPS::std_basis<3>(0),
                                DROPS::std_basis<3>(1),
                                DROPS::std_basis<3>(2),
                                DROPS::std_basis<3>(3),
                                20, 20, 20);
    MultiGridCL mg( brick);

    LevelsetP2CL lset( mg);

    lset.idx.CreateNumbering( 0, mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( &sphere_2);

    IdxDescCL ifaceidx( P1_FE);
    ifaceidx.SetTriangLevel( 0);
    IdxT num_unknowns = ifaceidx.NumUnknowns();
    CreateNumbOnInterfaceVertex( ifaceidx.GetIdx(), num_unknowns,
        ifaceidx.NumUnknownsVertex(), mg.GetTriangTetraBegin( 0), mg.GetTriangTetraEnd( 0),
        lset.Phi);
    ifaceidx.SetNumUnknowns( num_unknowns);
    std::cout << "NumUnknowns: " << ifaceidx.NumUnknowns() << std::endl;

    MatDescCL M( &ifaceidx, &ifaceidx);
    SetupInterfaceMassP1( mg, &M, lset.Phi);
    std::cerr << "Writing matrix to m_iface.txt\n";
    std::ofstream fff( "m_iface.txt");
    fff.precision( 18);
    fff << M.Data << std::endl;
    fff.close();    
   
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}
