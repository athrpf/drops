//**************************************************************************
// File:    levelset.tpp                                                   *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include <fstream>

namespace DROPS
{

inline double Sign( double x)
{
    return x<0 ? -1 : x>0 ? 1 : 0;
}

inline double SmoothedSign( double x, double alpha)
{
    return x/sqrt(x*x+alpha*alpha);
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::CreateNumbering(Uint level, IdxDescCL* idx)
// used for numbering of the Unknowns depending on the index IdxDesc[idxnum].
// sets up the description of the index idxnum in IdxDesc[idxnum],
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, that IdxDesc[idxnum].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx->TriangLevel = level;
    idx->NumUnknowns = 0;

    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL

    // allocate space for indices; number unknowns in TriangLevel level
    CreateNumbOnVertex( idxnum, idx->NumUnknowns, idx->NumUnknownsVertex,
                        _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level),
                        _dummyBnd);
    CreateNumbOnEdge( idxnum, idx->NumUnknowns, idx->NumUnknownsEdge,
                      _MG.GetTriangEdgeBegin(level), _MG.GetTriangEdgeEnd(level),
                      _dummyBnd );
}


template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::DeleteNumbering( IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetAllVertexBegin(level), _MG.GetAllVertexEnd(level) );
    DeleteNumbOnSimplex( idxnum, _MG.GetAllEdgeBegin(level), _MG.GetAllEdgeEnd(level) );
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::Init( scalar_fun_ptr phi0)
{
    const Uint lvl= Phi.RowIdx->TriangLevel,
               idx= Phi.RowIdx->GetIdx();
	 
    for (MultiGridCL::TriangVertexIteratorCL it= _MG.GetTriangVertexBegin(lvl),
        end= _MG.GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( it->GetCoord());
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= _MG.GetTriangEdgeBegin(lvl),
        end= _MG.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( GetBaryCenter( *it));
    }
}

inline double QuadVel( double f[10], int i)
{
    double sum, result;
    if (i<4) // hat function on vert
    {
        // Q = sum c[i]*f[i]
        // Gewichte c[i] = 1/420 	fuer vert i
        //                 1/2520	fuer uebrige verts
        //                -1/630        fuer an vert i anliegende edges
        //                -1/420        fuer uebrige drei edges   
        result= f[i]*(1/420.-1./2520.);
        sum= 0.;
        for (int k=0; k<4; ++k)
            sum+= f[k];
        result+= sum/2520.;
        sum= 0.;
        for (int k=0; k<3; ++k)
            sum+= f[EdgeByVert(i,k)+4];
        result+= -sum/630.;
        sum= 0.;
        const int oppF=OppFace(i);
        for (int k=0; k<3; ++k)
            sum+= f[EdgeOfFace(oppF, k)+4];
        result+= -sum/420.;
        return result;
    }
    else  // hat function on edge
    {
        i-= 4;
        // Q = sum c[i]*f[i]
        // Gewichte c[i] = 4/315 	fuer egde i
        //                 1/315	fuer opposite edge
        //                 2/315	fuer uebrige edges
        //                -1/630    	fuer an edge i anliegende verts
        //                -1/420        fuer uebrige zwei verts   
        result=  f[i+4]*4./315.;
        const int opp= OppEdge(i);
        result+= f[opp+4]/315.;
        sum= 0.;
        for(int k=0; k<6; ++k)
            if (k!=i && k!=opp)
                sum+= f[k+4];
        result+= sum*2./315.;
        sum= f[VertOfEdge(i,0)] + f[VertOfEdge(i,1)];
        result+= -sum/630.;
        sum= f[VertOfEdge(OppEdge(i),0)] + f[VertOfEdge(OppEdge(i),1)];
        result+= -sum/420.;

        return result;
    }
}

inline double GetMassP2( int i, int j)
{
// zur Erlaeuterung der zurueckgegebenen Gewichte
// beachte Kommentare zu Funktion QuadVel obendrueber!

    if (i>j) { int h=i; i=j; j=h; } // swap such that i<=j holds
    if (j<4) // i,j are vertex-indices
    {
        if (i==j)
            return 1./420.;
        else
            return 1./2520.;
    }
    else if (i>=4) // i,j are edge-indices
    {
        if (i==j)
            return 4./315.;
        else
        {
            if (i==OppEdge(j))
                return 1./315.;
            else
                return 2./315.;
        }
    }
    else // i<4, j>=4: Sonderfall...
    {
        if (i==VertOfEdge(j,0) || i==VertOfEdge(j,1))
            return -1./630.;
        else
            return -1./420.;
    }
}


inline double QuadVelGrad( double f[10], double g[4], int phi_i)
// approximates int ( f*g*v_i ) on reference tetrahedron,
// where f, g are scalar functions. 
// f is P2, g is P1 and v_i is the P2 ansatz function.
{
    switch( phi_i)
    {
      case 0: return -f[8] * g[3] / 0.1260e4 - f[7] * g[3] / 0.1260e4 - f[9] * g[0] / 0.2520e4 - f[9] * g[1] / 0.2520e4 - f[6] * g[2] / 0.1260e4 - f[5] * g[3] / 0.2520e4 - f[5] * g[2] / 0.1260e4 + f[2] * g[3] / 0.5040e4 - f[6] * g[3] / 0.2520e4 + f[3] * g[1] / 0.5040e4 + f[3] * g[2] / 0.5040e4 - f[7] * g[2] / 0.2520e4 - f[8] * g[2] / 0.2520e4 - f[4] * g[3] / 0.2520e4 - f[4] * g[2] / 0.2520e4 - f[9] * g[3] / 0.1260e4 + f[0] * g[3] / 0.2520e4 + f[1] * g[3] / 0.5040e4 + f[0] * g[2] / 0.2520e4 + f[1] * g[2] / 0.5040e4 - f[7] * g[1] / 0.2520e4 - f[9] * g[2] / 0.1260e4 - g[0] * f[8] / 0.2520e4 - f[8] * g[1] / 0.1260e4 - g[0] * f[6] / 0.2520e4 - f[6] * g[1] / 0.1260e4 + f[2] * g[1] / 0.5040e4 - f[5] * g[1] / 0.2520e4 + f[0] * g[1] / 0.2520e4 - f[4] * g[1] / 0.1260e4 + f[0] * g[0] / 0.840e3;
      case 1: return -f[8] * g[3] / 0.1260e4 - f[7] * g[3] / 0.1260e4 - f[9] * g[0] / 0.2520e4 - f[9] * g[1] / 0.2520e4 - f[6] * g[2] / 0.1260e4 - f[5] * g[3] / 0.2520e4 - f[5] * g[2] / 0.1260e4 + f[2] * g[3] / 0.5040e4 - f[6] * g[3] / 0.2520e4 + f[3] * g[2] / 0.5040e4 - f[7] * g[2] / 0.2520e4 + g[0] * f[3] / 0.5040e4 - f[8] * g[2] / 0.2520e4 - f[4] * g[3] / 0.2520e4 - f[4] * g[2] / 0.2520e4 - f[9] * g[3] / 0.1260e4 + f[0] * g[3] / 0.5040e4 + f[1] * g[3] / 0.2520e4 + f[0] * g[2] / 0.5040e4 + f[1] * g[2] / 0.2520e4 - f[7] * g[1] / 0.2520e4 - f[9] * g[2] / 0.1260e4 - g[0] * f[8] / 0.2520e4 - g[0] * f[7] / 0.1260e4 - g[0] * f[6] / 0.2520e4 - f[5] * g[1] / 0.2520e4 + g[0] * f[2] / 0.5040e4 - g[0] * f[5] / 0.1260e4 + g[0] * f[1] / 0.2520e4 - g[0] * f[4] / 0.1260e4 + f[1] * g[1] / 0.840e3;
      case 2: return -f[6] * g[1] / 0.1260e4 - f[8] * g[3] / 0.1260e4 - f[7] * g[3] / 0.1260e4 - f[9] * g[0] / 0.2520e4 - f[9] * g[1] / 0.2520e4 + f[2] * g[2] / 0.840e3 - f[5] * g[3] / 0.2520e4 + f[2] * g[3] / 0.2520e4 - f[6] * g[3] / 0.2520e4 + f[3] * g[1] / 0.5040e4 - f[7] * g[2] / 0.2520e4 + g[0] * f[3] / 0.5040e4 - f[8] * g[2] / 0.2520e4 - f[4] * g[3] / 0.2520e4 - f[4] * g[2] / 0.2520e4 - f[9] * g[3] / 0.1260e4 + f[0] * g[3] / 0.5040e4 + f[1] * g[3] / 0.5040e4 - f[7] * g[1] / 0.2520e4 - g[0] * f[8] / 0.2520e4 - g[0] * f[7] / 0.1260e4 - g[0] * f[6] / 0.2520e4 + f[2] * g[1] / 0.2520e4 - f[5] * g[1] / 0.2520e4 + g[0] * f[2] / 0.2520e4 - f[8] * g[1] / 0.1260e4 - g[0] * f[5] / 0.1260e4 + f[0] * g[1] / 0.5040e4 + g[0] * f[1] / 0.5040e4 - g[0] * f[4] / 0.1260e4 - f[4] * g[1] / 0.1260e4;
      case 3: return -f[6] * g[1] / 0.1260e4 + f[3] * g[3] / 0.840e3 - f[9] * g[0] / 0.2520e4 - f[9] * g[1] / 0.2520e4 - f[6] * g[2] / 0.1260e4 - f[5] * g[3] / 0.2520e4 - f[5] * g[2] / 0.1260e4 - f[6] * g[3] / 0.2520e4 + f[3] * g[1] / 0.2520e4 + f[3] * g[2] / 0.2520e4 - f[7] * g[2] / 0.2520e4 + g[0] * f[3] / 0.2520e4 - f[8] * g[2] / 0.2520e4 - f[4] * g[3] / 0.2520e4 - f[4] * g[2] / 0.2520e4 + f[0] * g[2] / 0.5040e4 + f[1] * g[2] / 0.5040e4 - f[7] * g[1] / 0.2520e4 - f[9] * g[2] / 0.1260e4 - g[0] * f[8] / 0.2520e4 - g[0] * f[7] / 0.1260e4 - g[0] * f[6] / 0.2520e4 + f[2] * g[1] / 0.5040e4 - f[5] * g[1] / 0.2520e4 + g[0] * f[2] / 0.5040e4 - f[8] * g[1] / 0.1260e4 - g[0] * f[5] / 0.1260e4 + f[0] * g[1] / 0.5040e4 + g[0] * f[1] / 0.5040e4 - g[0] * f[4] / 0.1260e4 - f[4] * g[1] / 0.1260e4;
      case 4: return f[6] * g[1] / 0.420e3 + f[7] * g[3] / 0.630e3 - f[3] * g[3] / 0.2520e4 + f[9] * g[0] / 0.1260e4 + f[9] * g[1] / 0.1260e4 + f[6] * g[2] / 0.630e3 - f[2] * g[2] / 0.2520e4 + f[5] * g[3] / 0.1260e4 + f[5] * g[2] / 0.630e3 - f[2] * g[3] / 0.2520e4 + f[6] * g[3] / 0.1260e4 - f[3] * g[1] / 0.1260e4 - f[3] * g[2] / 0.2520e4 + f[7] * g[2] / 0.1260e4 - g[0] * f[3] / 0.1260e4 + f[8] * g[2] / 0.1260e4 + f[4] * g[3] / 0.630e3 + f[4] * g[2] / 0.630e3 - f[0] * g[3] / 0.2520e4 - f[1] * g[3] / 0.2520e4 - f[0] * g[2] / 0.2520e4 - f[1] * g[2] / 0.2520e4 + f[7] * g[1] / 0.630e3 + f[9] * g[2] / 0.1260e4 + g[0] * f[8] / 0.630e3 + g[0] * f[7] / 0.420e3 + g[0] * f[6] / 0.630e3 - f[2] * g[1] / 0.1260e4 + f[5] * g[1] / 0.630e3 - g[0] * f[2] / 0.1260e4 + f[8] * g[1] / 0.420e3 + g[0] * f[5] / 0.420e3 - f[0] * g[1] / 0.1260e4 - g[0] * f[1] / 0.1260e4 + g[0] * f[4] / 0.210e3 + f[4] * g[1] / 0.210e3 + f[8] * g[3] / 0.630e3 + f[9] * g[3] / 0.1260e4;
      case 5: return f[6] * g[1] / 0.630e3 + f[7] * g[3] / 0.630e3 - f[3] * g[3] / 0.2520e4 + f[9] * g[0] / 0.630e3 + f[9] * g[1] / 0.1260e4 + f[6] * g[2] / 0.420e3 + f[5] * g[3] / 0.630e3 + f[5] * g[2] / 0.210e3 - f[2] * g[3] / 0.2520e4 + f[6] * g[3] / 0.1260e4 - f[3] * g[1] / 0.2520e4 - f[3] * g[2] / 0.1260e4 + f[7] * g[2] / 0.630e3 - g[0] * f[3] / 0.1260e4 + f[8] * g[2] / 0.1260e4 + f[4] * g[3] / 0.1260e4 + f[4] * g[2] / 0.630e3 - f[0] * g[3] / 0.2520e4 - f[1] * g[3] / 0.2520e4 - f[0] * g[2] / 0.1260e4 - f[1] * g[2] / 0.1260e4 + f[7] * g[1] / 0.1260e4 + f[9] * g[2] / 0.420e3 + g[0] * f[8] / 0.1260e4 + g[0] * f[7] / 0.420e3 + g[0] * f[6] / 0.630e3 - f[2] * g[1] / 0.2520e4 + f[5] * g[1] / 0.630e3 - g[0] * f[2] / 0.1260e4 + f[8] * g[1] / 0.1260e4 + g[0] * f[5] / 0.210e3 - f[0] * g[1] / 0.2520e4 - g[0] * f[1] / 0.1260e4 + g[0] * f[4] / 0.420e3 - f[1] * g[1] / 0.2520e4 + f[4] * g[1] / 0.630e3 + f[8] * g[3] / 0.1260e4 + f[9] * g[3] / 0.630e3;
      case 6: return f[6] * g[1] / 0.210e3 + f[7] * g[3] / 0.1260e4 - f[3] * g[3] / 0.2520e4 + f[9] * g[0] / 0.1260e4 + f[9] * g[1] / 0.630e3 + f[6] * g[2] / 0.210e3 + f[5] * g[3] / 0.1260e4 + f[5] * g[2] / 0.420e3 - f[2] * g[3] / 0.2520e4 + f[6] * g[3] / 0.630e3 - f[3] * g[1] / 0.1260e4 - f[3] * g[2] / 0.1260e4 + f[7] * g[2] / 0.1260e4 - g[0] * f[3] / 0.2520e4 + f[8] * g[2] / 0.630e3 + f[4] * g[3] / 0.1260e4 + f[4] * g[2] / 0.630e3 - f[0] * g[3] / 0.2520e4 - f[1] * g[3] / 0.2520e4 - f[0] * g[2] / 0.1260e4 - f[1] * g[2] / 0.1260e4 + f[7] * g[1] / 0.1260e4 + f[9] * g[2] / 0.420e3 + g[0] * f[8] / 0.1260e4 + g[0] * f[7] / 0.1260e4 + g[0] * f[6] / 0.630e3 - f[2] * g[1] / 0.1260e4 + f[5] * g[1] / 0.630e3 - g[0] * f[2] / 0.2520e4 + f[8] * g[1] / 0.420e3 + g[0] * f[5] / 0.630e3 - f[0] * g[1] / 0.1260e4 - g[0] * f[1] / 0.2520e4 + g[0] * f[4] / 0.630e3 + f[4] * g[1] / 0.420e3 + f[8] * g[3] / 0.630e3 + f[9] * g[3] / 0.630e3 - f[0] * g[0] / 0.2520e4;
      case 7: return f[6] * g[1] / 0.1260e4 + f[7] * g[3] / 0.210e3 + f[9] * g[0] / 0.630e3 + f[9] * g[1] / 0.1260e4 + f[6] * g[2] / 0.1260e4 - f[2] * g[2] / 0.2520e4 + f[5] * g[3] / 0.630e3 + f[5] * g[2] / 0.630e3 - f[2] * g[3] / 0.1260e4 + f[6] * g[3] / 0.1260e4 - f[3] * g[1] / 0.2520e4 - f[3] * g[2] / 0.2520e4 + f[7] * g[2] / 0.630e3 - g[0] * f[3] / 0.1260e4 + f[8] * g[2] / 0.1260e4 + f[4] * g[3] / 0.630e3 + f[4] * g[2] / 0.1260e4 - f[0] * g[3] / 0.1260e4 - f[1] * g[3] / 0.1260e4 - f[0] * g[2] / 0.2520e4 - f[1] * g[2] / 0.2520e4 + f[7] * g[1] / 0.630e3 + f[9] * g[2] / 0.630e3 + g[0] * f[8] / 0.630e3 + g[0] * f[7] / 0.210e3 + g[0] * f[6] / 0.1260e4 - f[2] * g[1] / 0.2520e4 + f[5] * g[1] / 0.1260e4 - g[0] * f[2] / 0.1260e4 + f[8] * g[1] / 0.630e3 + g[0] * f[5] / 0.420e3 - f[0] * g[1] / 0.2520e4 - g[0] * f[1] / 0.1260e4 + g[0] * f[4] / 0.420e3 - f[1] * g[1] / 0.2520e4 + f[4] * g[1] / 0.630e3 + f[8] * g[3] / 0.420e3 + f[9] * g[3] / 0.420e3;
      case 8: return f[6] * g[1] / 0.420e3 + f[7] * g[3] / 0.420e3 + f[9] * g[0] / 0.1260e4 + f[9] * g[1] / 0.630e3 + f[6] * g[2] / 0.630e3 - f[2] * g[2] / 0.2520e4 + f[5] * g[3] / 0.1260e4 + f[5] * g[2] / 0.1260e4 - f[2] * g[3] / 0.1260e4 + f[6] * g[3] / 0.630e3 - f[3] * g[1] / 0.1260e4 - f[3] * g[2] / 0.2520e4 + f[7] * g[2] / 0.1260e4 - g[0] * f[3] / 0.2520e4 + f[8] * g[2] / 0.630e3 + f[4] * g[3] / 0.630e3 + f[4] * g[2] / 0.1260e4 - f[0] * g[3] / 0.1260e4 - f[1] * g[3] / 0.1260e4 - f[0] * g[2] / 0.2520e4 - f[1] * g[2] / 0.2520e4 + f[7] * g[1] / 0.630e3 + f[9] * g[2] / 0.630e3 + g[0] * f[8] / 0.630e3 + g[0] * f[7] / 0.630e3 + g[0] * f[6] / 0.1260e4 - f[2] * g[1] / 0.1260e4 + f[5] * g[1] / 0.1260e4 - g[0] * f[2] / 0.2520e4 + f[8] * g[1] / 0.210e3 + g[0] * f[5] / 0.1260e4 - f[0] * g[1] / 0.1260e4 - g[0] * f[1] / 0.2520e4 + g[0] * f[4] / 0.630e3 + f[4] * g[1] / 0.420e3 + f[8] * g[3] / 0.210e3 + f[9] * g[3] / 0.420e3 - f[0] * g[0] / 0.2520e4;
      case 9: return f[6] * g[1] / 0.630e3 + f[7] * g[3] / 0.420e3 + f[9] * g[0] / 0.630e3 + f[9] * g[1] / 0.630e3 + f[6] * g[2] / 0.420e3 + f[5] * g[3] / 0.630e3 + f[5] * g[2] / 0.420e3 - f[2] * g[3] / 0.1260e4 + f[6] * g[3] / 0.630e3 - f[3] * g[1] / 0.2520e4 - f[3] * g[2] / 0.1260e4 + f[7] * g[2] / 0.630e3 - g[0] * f[3] / 0.2520e4 + f[8] * g[2] / 0.630e3 + f[4] * g[3] / 0.1260e4 + f[4] * g[2] / 0.1260e4 - f[0] * g[3] / 0.1260e4 - f[1] * g[3] / 0.1260e4 - f[0] * g[2] / 0.1260e4 - f[1] * g[2] / 0.1260e4 + f[7] * g[1] / 0.1260e4 + f[9] * g[2] / 0.210e3 + g[0] * f[8] / 0.1260e4 + g[0] * f[7] / 0.630e3 + g[0] * f[6] / 0.1260e4 - f[2] * g[1] / 0.2520e4 + f[5] * g[1] / 0.1260e4 - g[0] * f[2] / 0.2520e4 + f[8] * g[1] / 0.630e3 + g[0] * f[5] / 0.630e3 - f[0] * g[1] / 0.2520e4 - g[0] * f[1] / 0.2520e4 + g[0] * f[4] / 0.1260e4 + f[4] * g[1] / 0.1260e4 + f[8] * g[3] / 0.420e3 + f[9] * g[3] / 0.210e3 - f[1] * g[1] / 0.2520e4 - f[0] * g[0] / 0.2520e4;
    }
    return 0;
}


inline double QuadVelGrad( SVectorCL<3> f[10], SMatrixCL<3,5> g, int phi_i)
// approximates int ( f*g*v_i ) on reference tetrahedron,
// where f, g are vector-valued functions.
// f is P2, g is P1 and v_i is the P2 ansatz function.
{
    double ff[10], gg[4], sum= 0;
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<10; ++j) ff[j]= f[j][i];
        for (int j=0; j<4; ++j)  gg[j]= g(i,j);
        sum+= QuadVelGrad( ff, gg, phi_i);
    }
    return sum;
}


inline double QuadVelGrad( SVectorCL<3> f[10], SMatrixCL<3,5> g, double h[10])
// approximates int ( f*g h ) on reference tetrahedron,
// where f, g are vector-valued functions, h is a scalar function. 
// f, h are P2, g is P1.
{
    double sum= 0;
    for (int i=0; i<10; ++i)
        sum+= h[i]* QuadVelGrad( f, g, i);
    return sum;
}


inline double QuadVelGrad( SVectorCL<3> f[10], SMatrixCL<3,5> g, SMatrixCL<3,5> h)
// approximates int ( f*g f*h ) on reference tetrahedron,
// where f, g, h are vector-valued functions.
// f is P2, g, h are P1.
{
    double sum= 0;
    
    for (int i=0; i<10; ++i)
    {
        SVectorCL<3> grad;
        if (i<4)
        {
            for (int k=0; k<3; ++k) 
                grad[k]= h(k,i);
        }
        else
        {
            for (int k=0; k<3; ++k) 
                grad[k]= 0.5*( h(k,VertOfEdge(i-4,0)) + h(k,VertOfEdge(i-4,1)) );
        }
        double fh_i= 0;
        for (int k=0; k<3; ++k)
            fh_i+= grad[k]*f[i][k];
        // fh_i = f*h in DoF i
        sum+= fh_i * QuadVelGrad( f, g, i);
    }
    
    return sum;
}
/*
template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::SetupSystem( const DiscVelSolCL& vel)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.RowIdx->TriangLevel,
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> E(&_E, num_unks, num_unks), 
                               H(&_H, num_unks, num_unks);
    IdxT Numb[10];
    
    std::cerr << "entering SetupSystem: " << num_unks << " levelset unknowns." << std::endl;                            

    // fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10]; // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    SVectorCL<3> u_loc[10];
    double det, absdet, h_T;

    GetGradientsOnRef(GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        MakeGradients( Grad, GradRef, T);
        absdet= fabs( det);
        h_T= pow( absdet, 1./3);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and u_loc
        for(int i=0; i<4; ++i)
        {
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
	    u_loc[i]= vel.val( *sit->GetVertex(i));
        }
        for(int i=0; i<6; ++i)
        {
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
	    u_loc[i+4]= vel.val( *sit->GetEdge(i));
        }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            for(int j=0; j<10; ++j)
            {
                E( Numb[i], Numb[j])+= GetMassP2(i,j)*absdet; 
                E( Numb[i], Numb[j])+= _SD*h_T*QuadVelGrad(u_loc,Grad[i],j)*absdet; 
                
                H( Numb[i], Numb[j])+= QuadVelGrad(u_loc,Grad[j],i)*absdet;
                H( Numb[i], Numb[j])+= _SD*h_T*QuadVelGrad(u_loc,Grad[j],Grad[i])*absdet;
            }
    }
    
    E.Build();
    H.Build();
    std::cerr << _E.num_nonzeros() << " nonzeros in E, "
              << _H.num_nonzeros() << " nonzeros in H! " << std::endl;
}
*/
template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::SetupSystem( const DiscVelSolCL& vel)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.RowIdx->TriangLevel,
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> E(&_E, num_unks, num_unks), 
                               H(&_H, num_unks, num_unks);
    IdxT Numb[10];
    
    std::cerr << "entering SetupSystem: " << num_unks << " levelset unknowns." << std::endl;                            

    // fill value part of matrices
    Quad2CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad2CL<double> u_Grad[10]; // fuer u grad v_i
    SMatrixCL<3,3> T;
    double det, absdet, h_T;

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        h_T= pow( absdet, 1./3.);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and u_loc
        for(int i=0; i<4; ++i)
        {
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
	    u_loc.val[i]= vel.val( *sit->GetVertex(i));
        }
        for(int i=0; i<6; ++i)
        {
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }
        u_loc.val[4]= vel.val( *sit, 0.25, 0.25, 0.25);

        for(int i=0; i<10; ++i)
            u_Grad[i]= dot( u_loc, Grad[i]);
        
        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            for(int j=0; j<10; ++j)
            {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
                E( Numb[i], Numb[j])+= ( GetMassP2(i,j) 
                                       + u_Grad[i].quadP2(j)*_SD*h_T )*absdet; 
                
                // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                H( Numb[i], Numb[j])+= ( u_Grad[j].quadP2(i)
                                       + Quad2CL<>(u_Grad[i]*u_Grad[j]).quad() * _SD*h_T )*absdet;
            }
    }
    
    E.Build();
    H.Build();
    std::cerr << _E.num_nonzeros() << " nonzeros in E, "
              << _H.num_nonzeros() << " nonzeros in H! " << std::endl;
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::Reparam( Uint steps, double dt)
// Reparametrization of the levelset function Phi
{
    VectorCL Psi= Phi.Data, b;
    MatrixCL L, R;

    for (Uint i=0; i<steps; ++i)
    {
        SetupReparamSystem( R, Psi, b);
        L.LinComb( 1., _E, dt*_theta, R);
        
        b*= dt;
        b+= _E*Psi - dt*(1.-_theta) * (R*Psi);
        _gm.Solve( L, Psi, b);
        std::cout << "Reparam: res = " << _gm.GetResid() << ", iter = " << _gm.GetIter() << std::endl;
    }
    
    Phi.Data= Psi;
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::SetupReparamSystem( MatrixCL& _R, const VectorCL& Psi, VectorCL& b)
// R, b describe the following terms used for reparametrization:  
// b_i  = ( S(Phi0),         v_i + SD * w(Psi) grad v_i )
// R_ij = ( w(Psi) grad v_j, v_i + SD * w(Psi) grad v_i )
// where v_i, v_j denote the ansatz functions
// and w(Psi) = sign(Phi0) * grad Psi / |grad Psi| the scaled gradient of Psi
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.RowIdx->TriangLevel,
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> R(&_R, num_unks, num_unks);
    b.resize( 0);
    b.resize( num_unks);
    
    IdxT Numb[10];
    
    std::cerr << "entering SetupReparamSystem: " << num_unks << " levelset unknowns." << std::endl;                            

    // fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10]; // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    SVectorCL<3> w_loc[10], grad_Psi[4];
    double Sign_Phi[10];
    double det, absdet, h_T;
    const double alpha= 0.1;  // for smoothing of signum fct
    
    GetGradientsOnRef(GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        MakeGradients( Grad, GradRef, T);
        absdet= fabs( det);
        h_T= pow( absdet, 1./3);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb
        for (int i=0; i<4; ++i)
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        for (int i=0; i<6; ++i)
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        for (int i=0; i<4; ++i)
        {
            Sign_Phi[i]= SmoothedSign( Phi.Data[ Numb[i]], alpha);
            for (int k=0; k<3; ++k)
            {
                double sum= 0;
                for (int l=0; l<10; ++l)
                    sum+= Psi[ Numb[l]]*Grad[l](k,i);
                grad_Psi[i][k]= sum;
            }
	    w_loc[i]= Sign_Phi[i]*grad_Psi[i]/grad_Psi[i].norm();
        }
        for (int i=0; i<6; ++i)
        {
            Sign_Phi[i+4]= SmoothedSign( Phi.Data[ Numb[i+4]], alpha);
            // gradient on edge is the average of the gradient on the vertices, since gradient is linear
            SVectorCL<3> gr= 0.5*(grad_Psi[VertOfEdge(i,0)]+grad_Psi[VertOfEdge(i,1)]);
	    w_loc[i+4]= Sign_Phi[i+4]*gr/gr.norm();
        }

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            // b_i  = ( S(Phi0),         v_i + SD * w(Psi) grad v_i )
            b[ Numb[i]]+= QuadVel( Sign_Phi, i)*absdet;
            b[ Numb[i]]+= _SD*h_T*QuadVelGrad(w_loc,Grad[i], Sign_Phi)*absdet; 
            for(int j=0; j<10; ++j)
            {
                // R_ij = ( w(Psi) grad v_j, v_i + SD * w(Psi) grad v_i )
                R( Numb[i], Numb[j])+= QuadVelGrad(w_loc,Grad[j],i)*absdet;
                R( Numb[i], Numb[j])+= _SD*h_T*QuadVelGrad(w_loc,Grad[j],Grad[i])*absdet;
            }
        }
    }
    
    R.Build();
    std::cerr << _R.num_nonzeros() << " nonzeros in R!" << std::endl;
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::ComputeRhs( VectorCL& rhs) const
{
    rhs= _E*Phi.Data - _dt*(1-_theta) * (_H*Phi.Data);
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::DoStep( const VectorCL& rhs)
{
    _gm.Solve( _L, Phi.Data, rhs);
    std::cout << "res = " << _gm.GetResid() << ", iter = " << _gm.GetIter() <<std::endl;
}

template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::DoStep()
{
    VectorCL rhs( Phi.Data.size());
    ComputeRhs( rhs);
    DoStep( rhs);
}

template<class StokesProblemT>
bool LevelsetP2CL<StokesProblemT>::Intersects( const TetraCL& t) const
// Teste, ob Phi in allen DoFs das gleiche Vorzeichen hat
{
    const Uint idx= Phi.RowIdx->GetIdx();
    double PhiV0= Phi.Data[t.GetVertex(0)->Unknowns(idx)];
    
    for (int i=1; i<4; ++i)
        if( PhiV0*Phi.Data[t.GetVertex(i)->Unknowns(idx)] <= 0) return true;
    for (int i=0; i<6; ++i)
        if( PhiV0*Phi.Data[t.GetEdge(i)->Unknowns(idx)] <= 0) return true;

    return false;
}

/*
template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::Reparam2()
// Reparametrisierung, die den 0-Levelset fest laesst
// Idee: Phi auf T durch festes g=|grad Phi|>0 teilen -> 0-Levelset bleibt fest
{
    const size_t num_unks= Phi.RowIdx->NumUnknowns;
    const Uint   idx=      Phi.RowIdx->GetIdx();
    VectorBaseCL<bool> nearby( false, num_unks);  // is DoF near the levelset?
    VectorBaseCL<int>  num( 0, num_unks);         // how many tetras have this DoF?
    VectorCL           div( num_unks);            // sum of |grad Phi| on this tetras
    int                num_ls= 0;                 // how many tetras intersect the levelset?
    double             div_ls= 0;                 // sum of |grad Phi| on this tetras
    double      PhiLoc[10], det;
    int         sign[10];
    int         num_sign[3]; // - 0 + 
    IdxT        Numb[10];
    Point3DCL   grad;
    
    SMatrixCL<3,3> T;
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    for (MultiGridCL::TriangTetraIteratorCL it=_MG.GetTriangTetraBegin(), end=_MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            if (v<4) 
            {
                PhiLoc[v]= Phi.Data[it->GetVertex(v)->Unknowns(idx)];
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            }
            else
            {
                PhiLoc[v]= Phi.Data[it->GetEdge(v-4)->Unknowns(idx)];
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);
            }
            sign[v]= std::abs(PhiLoc[v])<1e-8 ? 0 : (PhiLoc[v]>0 ? 1 : -1);
        }
            
        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];
            const bool intersec= !(num_sign[0]*num_sign[2]==0 && num_sign[1]!=3); // (...) = no change of sign on child
            const double Phi0= PhiLoc[data.Vertices[0]];
            for (int i=0; i<3; ++i)
                grad[i]= PhiLoc[data.Vertices[i+1]] - Phi0;
            const double g= (T*grad).norm()*2; // doppelt so gross, da hier Kind des Tetra betrachtet wird

            if (intersec)
            {
                div_ls+= g;    ++num_ls;
                for (int vert= 0; vert<4; ++vert)
                    nearby[Numb[data.Vertices[vert]]]= true;
            }
            else
                for (int vert= 0; vert<4; ++vert)
                {
                    const size_t nv= Numb[data.Vertices[vert]];
                    div[nv]+= g;    ++num[nv];
                }
        }
    }
    
    const double c= num_ls/div_ls;
    std::cerr << "ReparamSaveIF: " << num_ls << " child tetras intersecting the levelset\n"
              << "\twith average norm of gradient = " << 1./c << std::endl;
    
    for (size_t i=0; i<num_unks; ++i)
        if (nearby[i])
            Phi.Data[i]*= c;
        else if (div[i]>1e-4)
        {
const double d= div[i]/num[i];
// Ausgabe nur bei rel. Abweichung von >50%
if (fabs(d-1./c)*c>0.5) std::cerr << d << "\t";
            Phi.Data[i]*= num[i]/div[i];
        }
std::cerr << std::endl;        
}
*/

inline void Solve2x2( const SMatrixCL<2,2>& A, SVectorCL<2>& x, const SVectorCL<2>& b)
{
    const double det= A(0,0)*A(1,1) - A(0,1)*A(1,0);
    x[0]= (A(1,1)*b[0]-A(0,1)*b[1])/det;
    x[1]= (A(0,0)*b[1]-A(1,0)*b[0])/det;
}


template<class StokesProblemT>
void LevelsetP2CL<StokesProblemT>::AccumulateBndIntegral( VecDescCL& f) const
{
    BaryCoordCL BaryDoF[10];
    Point3DCL   Coord[10];
    BaryDoF[0][0]= BaryDoF[1][1]= BaryDoF[2][2]= BaryDoF[3][3]= 1.;
    for (int edge=0; edge<6; ++edge)
        BaryDoF[edge+4]= 0.5*(BaryDoF[VertOfEdge(edge,0)] + BaryDoF[VertOfEdge(edge,1)]);

    const Uint  idx_phi= Phi.RowIdx->GetIdx(),
                idx_f=     f.RowIdx->GetIdx();
    double      PhiLoc[10];
    int         sign[10];
    int         num_sign[3]; // - 0 + 
    BaryCoordCL Bary[4];
    Point3DCL   PQRS[4], B[3];
    Point2DCL   AT_i, ab, tmp;
    IdxT        Numb[10];
/*    
std::ofstream fil("surf.off");
fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";
*/
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    P2DiscCL::GetGradientsOnRef( GradRef);
 
    for (MultiGridCL::TriangTetraIteratorCL it=_MG.GetTriangTetraBegin(), end=_MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, tmp[0], *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            Coord[v]= v<4 ? it->GetVertex(v)->GetCoord() : GetBaryCenter( *it->GetEdge(v-4));
            if (v<4) 
            {
                PhiLoc[v]= Phi.Data[it->GetVertex(v)->Unknowns(idx_phi)];
                Numb[v]= it->GetVertex(v)->Unknowns.Exist(idx_f) ?
                            it->GetVertex(v)->Unknowns(idx_f) : NoIdx;
            }
            else
            {
                PhiLoc[v]= Phi.Data[it->GetEdge(v-4)->Unknowns(idx_phi)];
                Numb[v]= it->GetEdge(v-4)->Unknowns.Exist(idx_f) ?
                             it->GetEdge(v-4)->Unknowns(idx_f) : NoIdx;
            }
            sign[v]= std::abs(PhiLoc[v])<1e-8 ? 0 : (PhiLoc[v]>0 ? 1 : -1);
        }
            
        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];
            if (num_sign[0]*num_sign[2]==0 && num_sign[1]!=3) // no change of sign on child
                continue;
            if (num_sign[1]==4)
            { 
                std::cerr << "WARNING: LevelsetP2CL::AccumulateBndIntegral: found 3-dim. zero level set, grid is too coarse!" << std::endl; 
                continue; 
            }

            int intersec= 0;
            // erst werden die Nullknoten in PQRS gespeichert...
            for (int vert= 0; vert<4; ++vert)
            {
                const int v= data.Vertices[vert];
                if (sign[v]==0)
                {
                    Bary[intersec]= BaryDoF[v];
                    PQRS[intersec++]= Coord[v];
                }
            }
            // ...dann die echten Schnittpunkte auf den Kanten mit Vorzeichenwechsel
            for (int edge= 0; edge<6; ++edge)
            {
                const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                          v1= data.Vertices[ VertOfEdge( edge, 1)];
                if (sign[v0]*sign[v1]<0) // different sign -> 0-level intersects this edge
                {
                    const double lambda= PhiLoc[v0]/(PhiLoc[v0]-PhiLoc[v1]);
                    Bary[intersec]= (1-lambda)*BaryDoF[v0] + lambda * BaryDoF[v1];
                    // bary-coords of tetra, not of subtetra!
                    PQRS[intersec++]= (1-lambda) * Coord[v0] + lambda * Coord[v1];
                }
            }
/*
fil << "geom {OFF " << intersec << " 1 0\n";
for (int i=0; i<intersec; ++i)
{
    for (int j=0; j<3; ++j)
        fil << PQRS[i][j] << ' ';
    fil << '\n';
}
if (intersec==4)
    fil << "4 0 1 3 2";
else
    fil << "3 0 1 2";
fil << "\n}\n";
*/

            if (intersec<3) continue; // Nullstellenmenge vom Mass 0!

            SMatrixCL<3,2> A;    // A = [ Q-P | R-P ]
            A(0,0)= PQRS[1][0]-PQRS[0][0];    A(0,1)= PQRS[2][0]-PQRS[0][0];
            A(1,0)= PQRS[1][1]-PQRS[0][1];    A(1,1)= PQRS[2][1]-PQRS[0][1];
            A(2,0)= PQRS[1][2]-PQRS[0][2];    A(2,1)= PQRS[2][2]-PQRS[0][2];
            SMatrixCL<2,2> ATA; 
            ATA(0,0)=           A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
            ATA(0,1)= ATA(1,0)= A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
            ATA(1,1)=           A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
            double sqrtDetATA= std::sqrt( ATA(0,0)*ATA(1,1) - ATA(1,0)*ATA(1,0) );
            BaryCoordCL BaryPQR, BarySQR;
            for (int i=0; i<3; ++i)
            {
                // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
                BaryPQR+= Bary[i];
                BarySQR+= Bary[i+1];
            
                // berechne B = A * (ATA)^-1 * AT
                AT_i[0]= A(i,0); AT_i[1]= A(i,1);
                Solve2x2(ATA,tmp,AT_i);
                B[i]= A*tmp;
            }

            if (intersec==4) // 4 intersections --> a+b != 1
            { // berechne a, b
                // Loese (Q-P)a + (R-P)b = S-P
                SMatrixCL<2,2> M;  
                M(0,0)= A(0,0); M(0,1)= A(0,1); M(1,0)=A(1,0); M(1,1)= A(1,1);
                // M is upper 2x2 part of A
                tmp[0]= PQRS[3][0]-PQRS[0][0]; tmp[1]= PQRS[3][1]-PQRS[0][1];
                // tmp = S-P
                Solve2x2( M, ab, tmp);

                //if (ab[0]<0 || ab[1]<0) 
                //    std::cerr<<"LevelsetP2CL::AccumulateBndIntegral: a or b negative"<<std::endl;
                // a,b>=0 muss erfuellt sein, da wegen edge+oppEdge==5 die Punkte P und S sich automatisch gegenueber liegen muessten...
            }

            sqrtDetATA/= 6;

            
            for (int v=0; v<10; ++v)
            {
                if (Numb[v]==NoIdx) continue;
                
                Point3DCL gr, gradv[4];
                // gr= grad Phi(P) + grad Phi(Q) + grad Phi(R)
                for (int node=0; node<4; ++node)
                { // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
                    gradv[node]= Grad[v].val[node];
                }
                    
                for (int k=0; k<4; ++k)
                    gr+= BaryPQR[k]*gradv[k];

                if (intersec==4)
                {
                    Point3DCL grSQR;
                    // grSQR = grad Phi(S) + grad Phi(Q) + grad Phi(R)
                    for (int k=0; k<4; ++k)
                        grSQR+= BarySQR[k]*gradv[k];

                    gr+= (ab[0] + ab[1] - 1) * grSQR;
                }
                // nun gilt:
                // gr = [grad Phi(P)+...] + (a+b-1)[grad Phi(S)+...]
                
                for (int i=0; i<3; ++i)
                {
                    const double val= inner_prod( gr, B[i]);
                    f.Data[Numb[v]+i]-= sigma * sqrtDetATA * val;
                }
            }
        } // Ende der for-Schleife ueber die Kinder
    }
//fil << "}\n";    
}

// ==============================================
//              CouplStokesLevelsetCL
// ==============================================

template <class StokesT, class SolverT>
CouplStokesLevelsetCL<StokesT,SolverT>::CouplStokesLevelsetCL
    ( StokesT& Stokes, LevelsetP2CL<StokesT>& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _rhs( Stokes.b.RowIdx->NumUnknowns), _ls_rhs( ls.Phi.RowIdx->NumUnknowns),
    _theta( theta)
{ 
    _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
    _old_curv->SetIdx( _b->RowIdx); _curv->SetIdx( _b->RowIdx);
    _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    ls.AccumulateBndIntegral( *_old_curv);
}

template <class StokesT, class SolverT>
CouplStokesLevelsetCL<StokesT,SolverT>::~CouplStokesLevelsetCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*(_old_b->Data + _old_curv->Data));
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _dt*_theta*_curv->Data, _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving Stokes took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    // setup system for levelset eq. as velocity has changed
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);
    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::CommitStep()
{
    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplStokesLevelsetCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_LvlSet.GetSolver().GetIter()==0) // no change of Phi -> no change of vel
            break;
    }
    CommitStep();
}


// ==============================================
//              CouplLevelsetStokesCL
// ==============================================

template <class StokesT, class SolverT>
CouplLevelsetStokesCL<StokesT,SolverT>::CouplLevelsetStokesCL
    ( StokesT& Stokes, LevelsetP2CL<StokesT>& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _rhs( Stokes.b.RowIdx->NumUnknowns), _ls_rhs( ls.Phi.RowIdx->NumUnknowns),
    _theta( theta)
{ 
    _old_b->SetIdx( _b->RowIdx); _cplM->SetIdx( _b->RowIdx); _old_cplM->SetIdx( _b->RowIdx);
    _old_curv->SetIdx( _b->RowIdx); _curv->SetIdx( _b->RowIdx);
    _Stokes.SetupInstatRhs( _old_b, &_Stokes.c, _old_cplM, _Stokes.t, _old_b, _Stokes.t);
    ls.AccumulateBndIntegral( *_old_curv);
}

template <class StokesT, class SolverT>
CouplLevelsetStokesCL<StokesT,SolverT>::~CouplLevelsetStokesCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupInstatRhs( _b, &_Stokes.c, _cplM, _Stokes.t, _b, _Stokes.t);

    _rhs=  _Stokes.A.Data * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data + _cplM->Data - _old_cplM->Data
         + _dt*( _theta*_b->Data + (1.-_theta)*(_old_b->Data + _old_curv->Data));
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    // setup system for levelset eq.
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);
    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _dt*_theta*_curv->Data, _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving Stokes took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::CommitStep()
{
    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplLevelsetStokesCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_solver.GetIter()==0) // no change of vel -> no change of Phi
            break;
    }
    CommitStep();
}


// ==============================================
//              CouplLevelsetStokes2PhaseCL
// ==============================================

template <class StokesT, class SolverT>
CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::CouplLevelsetStokes2PhaseCL
    ( StokesT& Stokes, LevelsetP2CL<StokesT>& ls, SolverT& solver, double theta)

  : _Stokes( Stokes), _solver( solver), _LvlSet( ls), _b( &Stokes.b), _old_b( new VelVecDescCL),
    _cplM( new VelVecDescCL), _old_cplM( new VelVecDescCL), 
    _curv( new VelVecDescCL), _old_curv( new VelVecDescCL), 
    _rhs( Stokes.v.RowIdx->NumUnknowns), _ls_rhs( ls.Phi.RowIdx->NumUnknowns),
    _theta( theta)
{ 
    Update(); 
}

template <class StokesT, class SolverT>
CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::~CouplLevelsetStokes2PhaseCL()
{
    if (_old_b == &_Stokes.b)
        delete _b;
    else
        delete _old_b; 
    delete _cplM; delete _old_cplM; delete _curv; delete _old_curv;
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::InitStep()
{
// compute all terms that don't change during the following FP iterations

    _LvlSet.ComputeRhs( _ls_rhs);

    _Stokes.t+= _dt;
    _Stokes.SetupRhs2( &_Stokes.c, _Stokes.t);

    _rhs=  _Stokes.A.Data * _Stokes.v.Data;
    _rhs*= (_theta-1.)*_dt;
    _rhs+= _Stokes.M.Data*_Stokes.v.Data - _old_cplM->Data
         + (_dt*(1.-_theta))*(_old_b->Data + _old_curv->Data);
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::DoFPIter()
// perform fixed point iteration
{
    TimerCL time;
    time.Reset();
    time.Start();

    // setup system for levelset eq.
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _LvlSet.SetTimeStep( _dt);

    time.Stop();
    std::cerr << "Discretizing Levelset took "<<time.GetTime()<<" sec.\n";
    time.Reset();

    _LvlSet.DoStep( _ls_rhs);

    time.Stop();
    std::cerr << "Solving Levelset took "<<time.GetTime()<<" sec.\n";

    time.Reset();
    time.Start();

    _curv->Clear();
    _LvlSet.AccumulateBndIntegral( *_curv);

    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _b, _cplM, _LvlSet, _Stokes.t);
    _mat.LinComb( 1., _Stokes.M.Data, _theta*_dt, _Stokes.A.Data);

    time.Stop();
    std::cerr << "Discretizing Stokes/Curv took "<<time.GetTime()<<" sec.\n";
    time.Reset();

    _Stokes.p.Data*= _dt;
    _solver.Solve( _mat, _Stokes.B.Data, _Stokes.v.Data, _Stokes.p.Data, 
                   _rhs + _cplM->Data + _dt*_theta*(_curv->Data + _b->Data), _Stokes.c.Data);
    _Stokes.p.Data/= _dt;

    time.Stop();
    std::cerr << "Solving Stokes took "<<time.GetTime()<<" sec.\n";
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::CommitStep()
{
    std::swap( _b, _old_b);
    std::swap( _cplM, _old_cplM);
    std::swap( _curv, _old_curv);
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::DoStep( int maxFPiter)
{
    if (maxFPiter==-1)
        maxFPiter= 999;

    InitStep();
    for (int i=0; i<maxFPiter; ++i)
    {
        DoFPIter();
        if (_solver.GetIter()==0 && _LvlSet.GetSolver().GetResid()<_LvlSet.GetSolver().GetTol()) // no change of vel -> no change of Phi
        {
            std::cerr << "Convergence after " << i+1 << " fixed point iterations!" << std::endl;
            break;
        }
    }
    CommitStep();
}

template <class StokesT, class SolverT>
void CouplLevelsetStokes2PhaseCL<StokesT,SolverT>::Update()
{
    IdxDescCL* const vidx= &_Stokes.vel_idx;
    IdxDescCL* const pidx= &_Stokes.pr_idx;
    TimerCL time;
    time.Reset();
    time.Start();

    std::cerr << "Updating discretization...\n";
    // IndexDesc setzen
    _b->SetIdx( vidx);       _old_b->SetIdx( vidx); 
    _cplM->SetIdx( vidx);    _old_cplM->SetIdx( vidx);
    _curv->SetIdx( vidx);    _old_curv->SetIdx( vidx);
    _rhs.resize( vidx->NumUnknowns);
    _ls_rhs.resize( _LvlSet.idx.NumUnknowns);
    _Stokes.c.SetIdx( pidx);
    _Stokes.A.SetIdx( vidx, vidx);
    _Stokes.B.SetIdx( pidx, vidx);
    _Stokes.M.SetIdx( vidx, vidx);

    // Diskretisierung
    _LvlSet.AccumulateBndIntegral( *_old_curv);
    _LvlSet.SetupSystem( _Stokes.GetVelSolution() );
    _Stokes.SetupSystem1( &_Stokes.A, &_Stokes.M, _old_b, _old_cplM, _LvlSet, _Stokes.t);
    _Stokes.SetupSystem2( &_Stokes.B, &_Stokes.c, _Stokes.t);
    
    time.Stop();
    std::cerr << "Discretizing took " << time.GetTime() << " sec.\n";
}

} // end of namespace DROPS

