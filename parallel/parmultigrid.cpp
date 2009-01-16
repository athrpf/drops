//****************************************************************************
// File:    parmultigrid.cpp                                                 *
// Content: Class that constitute the parallel multigrid                      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen      *
//          Oliver Fortmeier, RZ RWTH Aachen                                 *
// Version: 0.1                                                              *
// Date:                                                                     *
// Begin:   November, 14th, 2005                                             *
//****************************************************************************
/// \author Oliver Fortmeier
/// \file parmultigrid.cpp

#ifndef _DROPS_PARMULTIGRID_
#define _DROPS_PARMULTIGRID_

#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "num/solver.h"
#include <ddd.h>
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

namespace DROPS
{
/****************************************************************************
* I N I T I A L   S T A T I C   A S S I G N M E N T                         *
*****************************************************************************
* initial assignment of the static members of ParMultiGridCL                *
****************************************************************************/
DDD_IF  ParMultiGridCL::_EdgeIF =0;
DDD_IF  ParMultiGridCL::_TetraIF=0;
DDD_IF  ParMultiGridCL::_FaceIF =0;
DDD_IF ParMultiGridCL::NotMasterIF_=0;
DDD_IF ParMultiGridCL::_GhToMaTetra_IF=0;

DDD_TYPE ParMultiGridCL::_BndPtT         =0;        // DDD_TYPE "0" means invalide (see. ddd/typemgr)
DDD_TYPE ParMultiGridCL::_ChildPtrT      =0;

MultiGridCL*    ParMultiGridCL::_mg       =0;
MG_VertexContT* ParMultiGridCL::_VertCont =0;
MG_EdgeContT*   ParMultiGridCL::_EdgeCont =0;
MG_FaceContT*   ParMultiGridCL::_FaceCont =0;
MG_TetraContT*  ParMultiGridCL::_TetraCont=0;

ParMultiGridCL::TetraPCT ParMultiGridCL::ToHandleTetra_= TetraPCT();

ParMultiGridCL::VecDescPCT ParMultiGridCL::_VecDesc =   VecDescPCT();
ParMultiGridCL::BufferCT   ParMultiGridCL::_RecvBuf =   BufferCT();
ParMultiGridCL::ScalBndCT   ParMultiGridCL::_ScalBnd =  ScalBndCT();
ParMultiGridCL::VecBndCT    ParMultiGridCL::_VecBnd    = VecBndCT();

bool ParMultiGridCL::TransferMode   = false;
bool ParMultiGridCL::PrioChangeMode = false;
IdxT ParMultiGridCL::_RecvBufPos    = 0;

int ParMultiGridCL::_level           =-1;
bool ParMultiGridCL::_UnkOnSimplex[3]={false, false, false};

VecDescCL* ParMultiGridCL::_actualVec=0;

std::ostream *ParMultiGridCL::_os   = 0;
bool          ParMultiGridCL::_sane = true;

/****************************************************************************
* W R A P P E R                                                             *
*****************************************************************************
* Converts C++ functions to C functions, so DDD can use them correctly      *
****************************************************************************/
extern "C" void DeleteObjC(void* buffer, size_t size, int ddd_type) {ParMultiGridCL::DeleteObj(buffer,size,ddd_type);}

extern "C" void HandlerVDeleteC(DDD_OBJ obj) {ParMultiGridCL::HandlerDelete<VertexCL>(obj);}
extern "C" void HandlerEDeleteC(DDD_OBJ obj) {ParMultiGridCL::HandlerDelete<EdgeCL>(obj);}
extern "C" void HandlerFDeleteC(DDD_OBJ obj) {ParMultiGridCL::HandlerDelete<FaceCL>(obj);}
extern "C" void HandlerTDeleteC(DDD_OBJ obj) {ParMultiGridCL::HandlerDelete<TetraCL>(obj);}

extern "C" DDD_OBJ HandlerVConstructorC(size_t s, DDD_PRIO p, DDD_ATTR l) {return ParMultiGridCL::HandlerVConstructor(s,p,l);}
extern "C" DDD_OBJ HandlerEConstructorC(size_t s, DDD_PRIO p, DDD_ATTR l) {return ParMultiGridCL::HandlerEConstructor(s,p,l);}
extern "C" DDD_OBJ HandlerFConstructorC(size_t s, DDD_PRIO p, DDD_ATTR l) {return ParMultiGridCL::HandlerFConstructor(s,p,l);}
extern "C" DDD_OBJ HandlerTConstructorC(size_t s, DDD_PRIO p, DDD_ATTR l) {return ParMultiGridCL::HandlerTConstructor(s,p,l);}

extern "C" void HandlerVXferC( DDD_OBJ o, DDD_PROC to, DDD_PRIO p) {ParMultiGridCL::HandlerVXfer(o,to,p);}
extern "C" void HandlerEXferC( DDD_OBJ o, DDD_PROC to, DDD_PRIO p) {ParMultiGridCL::HandlerEXfer(o,to,p);}
extern "C" void HandlerFXferC( DDD_OBJ o, DDD_PROC to, DDD_PRIO p) {ParMultiGridCL::HandlerFXfer(o,to,p);}
extern "C" void HandlerTXferC( DDD_OBJ o, DDD_PROC to, DDD_PRIO p) {ParMultiGridCL::HandlerTXfer(o,to,p);}

extern "C" void HandlerVGatherC(DDD_OBJ o, int i, DDD_TYPE t, void* d) {ParMultiGridCL::HandlerVGather(o,i,t,d);}
extern "C" void HandlerEGatherC(DDD_OBJ o, int i, DDD_TYPE t, void* d) {ParMultiGridCL::HandlerEGather( o,i,t,d);}
extern "C" void HandlerTGatherC(DDD_OBJ o, int i, DDD_TYPE t, void* d) {ParMultiGridCL::HandlerTGather( o,i,t,d);}

extern "C" void HandlerVScatterC(DDD_OBJ o, int i, DDD_TYPE t, void* d, int n) {ParMultiGridCL::HandlerVScatter(o,i,t,d,n);}
extern "C" void HandlerEScatterC(DDD_OBJ o, int i, DDD_TYPE t, void* d, int n) {ParMultiGridCL::HandlerEScatter(o,i,t,d,n);}
extern "C" void HandlerTScatterC(DDD_OBJ o, int i, DDD_TYPE t, void* d, int n) {ParMultiGridCL::HandlerTScatter(o,i,t,d,n);}

extern "C" void HandlerTUpdateC(DDD_OBJ o)              {ParMultiGridCL::HandlerTUpdate(o);}
extern "C" void HandlerTObjMkConsC( DDD_OBJ o, int i)   {ParMultiGridCL::HandlerTObjMkCons(o,i);}
extern "C" void HandlerTSetPrioC(DDD_OBJ o, DDD_PRIO p) {ParMultiGridCL::HandlerTSetPrio(o,p);}

extern "C" int GatherEdgeMFRC(DDD_OBJ o, void* d) {return ParMultiGridCL::GatherEdgeMFR(o,d);}
extern "C" int ScatterEdgeMFRC(DDD_OBJ o, void* d) {return ParMultiGridCL::ScatterEdgeMFR(o,d);}

extern "C" int GatherTetraRestrictMarksC( DDD_OBJ o, void* b) {return ParMultiGridCL::GatherTetraRestrictMarks(o,b);}
extern "C" int ScatterTetraRestrictMarksC(DDD_OBJ o, void* b) {return ParMultiGridCL::ScatterTetraRestrictMarks(o,b);}

extern "C" int ExecGhostRescueC(         DDD_OBJ o) {return ParMultiGridCL::ExecGhostRescue(o);}
extern "C" int ExecGhVertRescueC(        DDD_OBJ o) {return ParMultiGridCL::ExecGhVertRescue(o);}
extern "C" int ExecHasGhostC(            DDD_OBJ o) {return ParMultiGridCL::ExecHasGhost(o);}
extern "C" int ExecAdaptVGhostMidVertexC(DDD_OBJ o) {return ParMultiGridCL::ExecAdaptVGhostMidVertex(o);}

extern "C" int GatherEdgeSaneC( DDD_OBJ o, void* d, DDD_PROC p, DDD_ATTR a) {return ParMultiGridCL::GatherEdgeSane(o,d,p,a);}
extern "C" int ScatterEdgeSaneC(DDD_OBJ o, void* d, DDD_PROC p, DDD_ATTR a) {return ParMultiGridCL::ScatterEdgeSane(o,d,p,a);}
extern "C" int GatherFaceSaneC( DDD_OBJ o, void* d, DDD_PROC p, DDD_ATTR a) {return ParMultiGridCL::GatherFaceSane(o,d,p,a);}
extern "C" int ScatterFaceSaneC(DDD_OBJ o, void* d, DDD_PROC p, DDD_ATTR a) {return ParMultiGridCL::ScatterFaceSane(o,d,p,a);}

extern "C" int GatherUnknownsRefC (DDD_OBJ o, void* b) { return ParMultiGridCL::GatherUnknownsRef(o,b); }
extern "C" int ScatterUnknownsRefC(DDD_OBJ o, void* b) { return ParMultiGridCL::ScatterUnknownsRef(o,b); }

extern "C" int GatherUnknownsMigVC (DDD_OBJ o, void* b) { return ParMultiGridCL::GatherUnknownsMigV(o,b); }
extern "C" int ScatterUnknownsMigVC(DDD_OBJ o, void* b) { return ParMultiGridCL::ScatterUnknownsMigV(o,b); }
extern "C" int GatherUnknownsMigEC (DDD_OBJ o, void* b) { return ParMultiGridCL::GatherUnknownsMigE(o,b); }
extern "C" int ScatterUnknownsMigEC(DDD_OBJ o, void* b) { return ParMultiGridCL::ScatterUnknownsMigE(o,b); }

extern "C" int GatherInterpolValuesVC (DDD_OBJ o, void* b) { return ParMultiGridCL::GatherInterpolValues<VertexCL>(o,b); }
extern "C" int ScatterInterpolValuesVC(DDD_OBJ o, void* b) { return ParMultiGridCL::ScatterInterpolValues<VertexCL>(o,b); }
extern "C" int GatherInterpolValuesEC (DDD_OBJ o, void* b) { return ParMultiGridCL::GatherInterpolValues<EdgeCL>(o,b); }
extern "C" int ScatterInterpolValuesEC(DDD_OBJ o, void* b) { return ParMultiGridCL::ScatterInterpolValues<EdgeCL>(o,b); }


/****************************************************************************
* C O N S T R U C T O R S                                                   *
*****************************************************************************
* Constructor an Destructor of the parallel multigrid                       *
* AttachTo: assign the given Multigrid to the parallel mutigrid             *
*           or tell the Multigrid about a VectorDescriber-Class             *
****************************************************************************/
/// \brief Constructor with number of Vector-Describer-Classes
///
/// Init the parallel stuff
ParMultiGridCL::ParMultiGridCL()
{
    Assert(EdgeCL::GetType()==0, DROPSErrCL("ParMultiGridCL: Constructor is called twice"),DebugParallelC);

    // Init the DDD stuff
    DeclareAll();       // Declare all DDD-Types
    DefineAll();        // Define all DDD-Types
    InitIF();           // Define all DDD-Interface
    SetAllHandler();    // Set the DDD-Handlers

    // Init the containers for the unknowns
    _RecvBuf.resize(0);

    // There are no active transfers neither Prio-Change-enviroment while creating this class
    TransferMode=false;
    PrioChangeMode=false;
    _UnkOnSimplex[0]=_UnkOnSimplex[1]=_UnkOnSimplex[2]=false;
}

/// \brief Constructor with a MGBuilderCL and the number of VecDescCL
///
/// see constructor above and assign the Multigrid given by the MGBuilderCL to this parallel MultiGrid
ParMultiGridCL::ParMultiGridCL(const MGBuilderCL& Builder)
{
    Assert(EdgeCL::GetType()==0, DROPSErrCL("ParMultiGridCL: Constructor is called twice"),DebugParallelC);

    // Init the DDD stuff
    DeclareAll();
    DefineAll();
    InitIF();
    SetAllHandler();

    // Create the multigrid
    _mg= new MultiGridCL( Builder);

    // assign the references to the simplex containers
    _VertCont=  &_mg->_Vertices;
    _EdgeCont=  &_mg->_Edges;
    _FaceCont=  &_mg->_Faces;
    _TetraCont= &_mg->_Tetras;

    // Init the containers for the unknowns
//     _SelectBnd.resize(numVecDesc, std::numeric_limits<size_t>::max());
    _RecvBuf.resize(0);

    // There are no active transfers neither Prio-Change-enviroment while creating this class
    TransferMode=false;
    PrioChangeMode=false;
    _UnkOnSimplex[0]=_UnkOnSimplex[1]=_UnkOnSimplex[2]=false;
}

ParMultiGridCL::~ParMultiGridCL()
{
//  if (_UnkOnSimplex)
//      delete[] _UnkOnSimplex;
}

/// \brief Assign the MultiGridCL to this ParMultiGridCL
void ParMultiGridCL::AttachTo( MultiGridCL& mg)
{
    // Store referenz to the Multigrid and to the simplex containers
    _mg= &mg;
    _VertCont=  &_mg->_Vertices;
    _EdgeCont=  &_mg->_Edges;
    _FaceCont=  &_mg->_Faces;
    _TetraCont= &_mg->_Tetras;
}

/// \brief Store a pointer to a scalar boundary condition, that belongs to an Index
template<>
  void ParMultiGridCL::AttachTo<BndDataCL<double> >(const IdxDescCL* idxDesc, const BndDataCL<double>* bndp)
/** Scalar boundary conditions are given. Hence, a pointer to the boundary
    consitions are stored _ScalBnd in at the same position where the VecDescCL
    can be found.
    \pre The VecDesc, corresponding to \a idxDesc,  must have been set by
         AttachTo(VecDescCL*) before calling this routine
    \param idxDesc IdxDescCL matching to \a bndp
    \param bndp    pointer to the correspondint BndDataCL  */
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert( vecPos!=_VecDesc.size() && vecPos<_ScalBnd.size(),
            DROPSErrCL("ParMultiGridCL::AttachTo<BndDataCL<double>>: VecDesc is not known so far, set it with AttachTo before calling this routine"),
            DebugParallelNumC);
    _ScalBnd[vecPos]= bndp;
}

/// \brief Store a pointer to a vectorial boundary condition, that belongs to an Index
template<>
  void ParMultiGridCL::AttachTo<BndDataCL<Point3DCL> >(const IdxDescCL* idxDesc, const BndDataCL<Point3DCL>* bndp)
/** Vectorial boundary conditions are given. Hence, a pointer to the boundary
    consitions are stored _VecBnd in at the same position where the VecDescCL
    can be found.
    \pre The VecDesc, corresponding to \a idxDesc,  must have been set by
         AttachTo(VecDescCL*) before calling this routine
    \param idxDesc IdxDescCL matching to \a bndp
    \param bndp    pointer to the correspondint BndDataCL  */
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert( vecPos!=_VecDesc.size() && vecPos<_VecBnd.size(),
            DROPSErrCL("ParMultiGridCL::AttachTo<BndDataCL<double>>: VecDesc is not known so far, set it with AttachTo before calling this routine"),
            DebugParallelNumC);
    _VecBnd[vecPos]= bndp;
}

/// \brief Delete all information about the VecDescCL
void ParMultiGridCL::DeleteVecDesc()
{
    _VecDesc.resize( 0);
    _VecBnd.resize( 0);
    _ScalBnd.resize( 0);
    _UnkOnSimplex[0]=false;
    _UnkOnSimplex[1]= false;
    _UnkOnSimplex[2]= false;
}

/// \brief Get scalar boundary condition to a known VecDesCL
template<>
  const BndDataCL<double>* ParMultiGridCL::GetBndCond<BndDataCL<double> >( const IdxDescCL* idxDesc)
/// First, find the position, where the VecDescCL is stored. Second, the pointer
/// to the boundary conditions is returned
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert(vecPos<_VecDesc.size() && vecPos<_ScalBnd.size() && _ScalBnd[vecPos]!=0,
           DROPSErrCL("ParMultiGridCL::GetBndCond<BndDataCL<double>>: BC not set"),
           DebugParallelNumC);
    return _ScalBnd[vecPos];
}

/// \brief Get vectorial boundary condition to a known VecDesCL
template<>
  const BndDataCL<Point3DCL>* ParMultiGridCL::GetBndCond<BndDataCL<Point3DCL> >( const IdxDescCL* idxDesc)
/// First, find the position, where the VecDescCL is stored. Second, the pointer
/// to the boundary conditions is returned
{
    size_t vecPos= GetStorePos( idxDesc);
    Assert(vecPos<_VecDesc.size() && vecPos<_VecBnd.size() && _VecBnd[vecPos]!=0,
           DROPSErrCL("ParMultiGridCL::GetBndCond<BndDataCL<Point3DCL>>: BC not set"),
           DebugParallelNumC);
    return _VecBnd[vecPos];
}

/// \brief Handle Unknowns after a refine operation
void ParMultiGridCL::HandleUnknownsAfterRefine(/*const std::vector<VecDescCL*> &newVecDesc*/)
/**
    This procedure performs the following operations:
    <ol>
     <li>Transfer unknowns from killed ghost tetras to the master tetras</li>
     <li>Delete all subsimplices, that are marked for removement by the refinement algorithm</li>
     <li>Unsubscribe the killed ghost tetras from the DDD-System</li>
     <li>Delete ghost tetras that are not needed any more</li>
     <li>Copy recieved values into the recieve buffer</li>
    </ol>
    \pre  The VecDescCL's, that describes the unknowns before the refinement algorithm has been performed,
          has to be attached to the ParMultiGridCL by the procedure 'AttachTo'.
    \post All subsimplices, that stores unknowns has an accumulated value of its unknowns
*/
{
    Assert(_RecvBuf.size()==0, DROPSErrCL("ParMultiGridCL::HandleUnknownsAfterRefine: Recieve Buffer is not empty!"), DebugParallelNumC);
    if (!_mg->UnknownsForRefine())
        return;
    Assert(VecDescRecv(), DROPSErrCL("ParMultiGridCL::HandleUnknownsAfterRefine: missing vector describers"), DebugParallelNumC);

    // Allocate mem for recieve buffer
    IdxT numKnownUnknowns=0;
    for (size_t i=0; i<_VecDesc.size(); ++i)
        numKnownUnknowns+=_VecDesc[i]->Data.size();
    _RecvBuf.resize(numKnownUnknowns);

    // Transfer the unknowns to the right processor
    DDD_IFOneway( _GhToMaTetra_IF,                               // transfer unknowns from killed ghost-tetra to corresponding master tetra
                  IF_FORWARD,                                    // transfer in one direction
                  NumberOfUnknownsOnTetra()*sizeof(TransferUnkT),// number of (double+bool) values
                  &GatherUnknownsRefC, &ScatterUnknownsRefC);    // handler for gather and scatter

    // Right now no more unused simplices are needed any more!
    _mg->PrepareModify();

    // kill all subsimplices, that has not beed done by refinement algorithm
    for (Uint l=0; l<=_mg->GetLastLevel(); ++l)
    {
        _mg->_Vertices[l].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
        _mg->_Edges[l].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
        _mg->_Faces[l].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    }

    // Kill Ghosts tetras that are not needed any more
    if (_mg->KilledGhosts())
    {
        std::list<MultiGridCL::TetraIterator>::iterator end(_mg->toDelGhosts_.end());
        // unsubscribe them from the DDD-System
        DDD_XferBegin();
        for (std::list<MultiGridCL::TetraIterator>::iterator it(_mg->toDelGhosts_.begin()); it!=end; ++it)
            (*it)->XferDelete();
        DDD_XferEnd();

        // physically delete them
        for (std::list<MultiGridCL::TetraIterator>::iterator it(_mg->toDelGhosts_.begin()); it!=end; ++it)
        {
            Assert((*it)->GetChildBegin()==0,
                   DROPSErrCL("ParMultiGridCL::HandleUnknownsAfterRefine: Killing Ghost tetra, that has children"),
                   DebugParallelNumC);
            _mg->_Tetras[(*it)->GetLevel()].erase(*it);
        }

        // information about ghost, that should be deleted are not needed any more
        _mg->toDelGhosts_.resize(0);
    }

    // there are no more ghost tetras for deleting
    _mg->killedGhostTetra_= false;

    // Adjust level
    while ( _mg->GetLastLevel()>0 && _mg->_Tetras[_mg->GetLastLevel()].empty() )
        _mg->RemoveLastLevel();
    AdjustLevel();

    // no more modifications has to be done
    _mg->FinalizeModify();
    _mg->ClearTriangCache();
}

/// \brief Fill a new vector describer class with datas (Has to be called after the refinement an migration algorithm!)
void ParMultiGridCL::HandleNewIdx(IdxDescCL* oldIdxDesc, VecDescCL* newVecDesc)
/** Within the refinement or migration algorithm a single processor
    can get new unknowns and give other unknowns to other processors. So each
    processor put the recieved unknowns in a buffer. In order to get a fully filled
    data vector, this routine builds this vector by taking values out of the
    recieve buffer or the "old" vector.*/
/*  This routine iterates over all vertices, edges and tetras and calls the routine PutData.
    PutData collects the data and put it into the new vector at the right position.*/
{
    Assert(VecDescRecv(), DROPSErrCL("ParMultiGridCL::HandleNewIdx: No Indices recieved before transfer"), DebugParallelNumC);

    // old and new index
    const Uint old_idx= oldIdxDesc->GetIdx(),
               new_idx= newVecDesc->RowIdx->GetIdx();

    // find the right index within _VecDesc
    size_t vecIdx= GetStorePos( oldIdxDesc);
    Assert(vecIdx!=_VecDesc.size(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: The appropriate VecDesc has not been set"),
           DebugParallelNumC);

    // create the new data vector
    newVecDesc->Data.resize(newVecDesc->RowIdx->NumUnknowns());

    // abbrevations for accessing the data vectors
    VectorCL *new_data= &(newVecDesc->Data);
    const VectorCL* const old_data= &(_VecDesc[vecIdx]->Data);

    // error checking
    Assert(newVecDesc->RowIdx->NumUnknownsVertex() == oldIdxDesc->NumUnknownsVertex(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: number of unknowns on vertices does not match"), DebugParallelNumC);
    Assert(newVecDesc->RowIdx->NumUnknownsEdge()==oldIdxDesc->NumUnknownsEdge(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: number of unknowns on edges does not match"), DebugParallelNumC);
    Assert(newVecDesc->RowIdx->NumUnknownsTetra()==oldIdxDesc->NumUnknownsTetra(),
           DROPSErrCL("ParMultiGridCL::HandleNewIdx: number of unknowns on tetras does not match"), DebugParallelNumC);

    // Put unknowns on vertices into the new vector
    for (MultiGridCL::const_VertexIterator sit=_mg->GetAllVertexBegin();
            sit!=_mg->GetAllVertexEnd(); ++sit)
    {
        PutData(sit, old_data, new_data, old_idx, new_idx, oldIdxDesc);
    }

    // Put unknowns on edges into the new vector and interpolate very special midvertices
    if (oldIdxDesc->NumUnknownsEdge()){
        for (MultiGridCL::const_EdgeIterator sit=_mg->GetAllEdgeBegin();
                sit!=_mg->GetAllEdgeEnd(); ++sit)
        {
            if (oldIdxDesc->NumUnknownsVertex()==1)
                PutData(sit, old_data, new_data, old_idx, new_idx, oldIdxDesc, GetBndCond<BndDataCL<double> >(oldIdxDesc));
            else
                PutData(sit, old_data, new_data, old_idx, new_idx, oldIdxDesc, GetBndCond<BndDataCL<Point3DCL> >(oldIdxDesc));
        }
    }
}

/// \brief Send interpolated values from processers that can interpolate a value to processors that could not interpolate a value
void ParMultiGridCL::CompleteRepair(VecDescCL* newVec)
/** After RepairAfterRefine[P1|P2] call this function to exchange interpolated
    values.
    \param newVec the interpolated VecDescCL*/
{

    Assert(newVec->RowIdx->NumUnknownsEdge()==0 || newVec->RowIdx->NumUnknownsEdge()==newVec->RowIdx->NumUnknownsVertex(),
           DROPSErrCL("ParMultiGridCL::CompleteRepair: Unknowns on vertices and edges must be the same"),
           DebugParallelNumC);

    if (newVec->RowIdx->NumUnknownsVertex()){
        _actualVec=newVec;             // gather and scatter functions has to know about the index and values
        DDD_IFExchange(InterfaceCL<VertexCL>::GetIF(),                              // exchange datas over distributed vertices
                    _actualVec->RowIdx->NumUnknownsVertex() * sizeof(TransferUnkT), // number of datas to be exchanged
                    GatherInterpolValuesVC,                                         // how to gather datas
                    ScatterInterpolValuesVC                                         // how to scatter datas
                    );
        _actualVec=0;                           // fortget about the actual index and values
    }

    if (newVec->RowIdx->NumUnknownsEdge()){
        _actualVec=newVec;
        DDD_IFExchange(InterfaceCL<EdgeCL>::GetIF(),
                    _actualVec->RowIdx->NumUnknownsEdge() * sizeof(TransferUnkT),
                    GatherInterpolValuesEC,
                    ScatterInterpolValuesEC
                    );
        _actualVec=0;
    }
}

/// \brief Put datas on a vertex from the old vector into a new vector
void ParMultiGridCL::PutData(MultiGridCL::const_VertexIterator& sit,
                             const VectorCL* const old_data, VectorCL* new_data,
                             const Uint old_idx, const Uint new_idx,
                             const IdxDescCL* idxDesc)
/** This routine puts the unknowns on a vertex according to an index into
    a new data vector. Therefore datas are taken from the recieve buffer, if the
    data has been stored on another proc before the refinement and migration, or
    out of the "old" vector, if the calling proc has already owned these data
    before the refinement.
    \param sit         Iterator onto a simplex
    \param old_data    Pointer to the old data vector
    \param new_data    Pointer to the new data vector
    \param old_idx     old index
    \param new_idx     new index
    \param idxDesc     describer of index, just used for getting information
                       about number of unknowns on simplices*/
{
    const Uint numUnknowns=idxDesc->NumUnknownsVertex();
    // check if there are unknowns to copy. (They have to exist on the new and old ParMultiGrid)
    if (sit->Unknowns.Exist()
        && sit->Unknowns.Exist(new_idx)
        && sit->Unknowns.Exist(old_idx) )
    {
        const IdxT new_sysnum = sit->Unknowns(new_idx),                 // position, where to find the new unknown(s)
                   old_sysnum = sit->Unknowns(old_idx);                 // position, where to find the old unknown(s)
        if (sit->Unknowns.UnkRecieved(old_idx))                         // this is an "new" unknown
        {
            for (Uint i=0; i<numUnknowns; ++i)
                (*new_data)[new_sysnum+i]= _RecvBuf[old_sysnum+i];
            sit->Unknowns.SetUnkRecieved(new_idx);                      // set recieved flag as a flag, that this unknown is an old one
        }
        else
        {
            for (Uint i=0; i<numUnknowns; ++i)
                (*new_data)[new_sysnum+i]=  (*old_data)[old_sysnum+i];
            sit->Unknowns.SetUnkRecieved(new_idx);                      // set recieved flag as a flag, that this unknown is an old one
        }
    }
}

/// \brief Save unknowns on edges, that are deleted within migration
void ParMultiGridCL::RescueUnknownsOnEdges()
/** Within the migration it may happen, that an edge is removed, but the
    midvertex stays on the processor. If the midvertex has unknowns after the
    migration on this processor and do not has unknowns before the refinement,
    the interpolation must be done, before the edge is deleted. So this is done
    here. <br>
    We distinguish between two cases:
    <ol>
     <li>P1-finite elements: The value on the midvertex is interpolated by the
         values on the two vertices of the edge. If a value is given by
         Dirichlet BC, this must be checked.</li>
     <li>P2-finite elements: The value on the edge is set to the midvertex. Here
         no check for Dirichlet BC is needed.</li>
    </ol>*/
{
    Assert(UnknownsOnSimplices(),
           DROPSErrCL("ParMultiGridCL::RescueUnknownsOnEdges: Called without unknowns"),
           DebugParallelNumC);

    // check all edges, if the edge is removed, but its midvertex stays on this
    // processor and the midvertex has no value
    for (int l=0; l<=_level; ++l){
        for (MultiGridCL::EdgeIterator sit((*_EdgeCont)[l].begin()), tEnd((*_EdgeCont)[l].end()); sit!=tEnd; ++sit){
            if (sit->IsRefined()
                && sit->IsMarkedForRemovement()
                && !sit->GetMidVertex()->IsMarkedForRemovement()
                && !sit->GetMidVertex()->Unknowns.Exist())
            {

                // rescue unknowns for all known indices
                for (size_t idx_type=0; idx_type<_VecDesc.size(); ++idx_type){
                    const VecDescCL* sol          = _VecDesc[idx_type];                 // Pointer to the VecDescCL
                    const Uint       idx          = sol->RowIdx->GetIdx();              // index of the unknowns
                    const Uint       numUnkOnEdge = sol->RowIdx->NumUnknownsEdge(),     // unknowns on edges
                                     numUnkOnVert = sol->RowIdx->NumUnknownsVertex();

                    if (numUnkOnEdge==0 && numUnkOnVert){                                 // Rescue unknowns for P1 Elements
                        if (numUnkOnVert==1){   // scalar unknowns
                            const BndDataCL<double>* bnd = GetBndCond<BndDataCL<double> >(sol->RowIdx);
                            double new_dof;
                            if(!bnd->IsOnDirBnd(*sit->GetMidVertex()) ){
                                if ( LinearInterpolation(*sit, idx, bnd, sol->Data, new_dof) ){
                                    PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx, new_dof);
                                    sit->GetMidVertex()->Unknowns.SetUnkRecieved(idx);
                                }
                                else
                                    throw DROPSErrCL("ParMultiGridCL::RescueUnknownsOnEdges: Cannot interpolate onto midvertex");
                            }
                        }
                        else{                   // vectorial unknowns
                            const BndDataCL<Point3DCL>* bnd = GetBndCond<BndDataCL<Point3DCL> >(sol->RowIdx);;
                            Point3DCL new_dof;
                            if(!bnd->IsOnDirBnd(*sit->GetMidVertex()) ){
                                if ( LinearInterpolation(*sit, idx, bnd, sol->Data, new_dof) ){
                                    PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx, new_dof);
                                    sit->GetMidVertex()->Unknowns.SetUnkRecieved(idx);
                                }
                                else
                                    throw DROPSErrCL("ParMultiGridCL::RescueUnknownsOnEdges: Cannot interpolate onto midvertex");
                            }
                        }
                    }
                    else{                                                               // Rescue unknowns for P2 Elements
                        if (sit->Unknowns.Exist()){
                            if (numUnkOnEdge==1)
                                PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx,
                                    GetDofOutOfVector<EdgeCL, VectorCL, double>()(*sit, idx, sol->Data));
                            else
                                PutDofIntoRecvBuffer(*sit->GetMidVertex(), idx,
                                    GetDofOutOfVector<EdgeCL, VectorCL, Point3DCL>()(*sit, idx, sol->Data));
                            sit->GetMidVertex()->Unknowns.SetUnkRecieved(idx);
                        }
                    }
                }
            }
        }
    }
}

/// \brief Delete recieve Buffer
/** After a refinement and a miragtion the recieve buffer is not needed any more */
void ParMultiGridCL::DeleteRecvBuffer()
{
    _RecvBuf.resize(0);
    _RecvBufPos=0;
}

/// \brief Get overall number of unknowns on a tetrahedron
/** Check number of unknowns on tetra and all its subsimplices for number of
    unknowns according to vecdesc */
Uint ParMultiGridCL::NumberOfUnknownsOnTetra()
{
    Uint result=0;
    for (VecDescPCT::const_iterator it(_VecDesc.begin()), end(_VecDesc.end()); it!=end; ++it)
        result +=  NumVertsC * (*it)->RowIdx->NumUnknownsVertex()
                 + NumEdgesC * (*it)->RowIdx->NumUnknownsEdge()
                 + 1         * (*it)->RowIdx->NumUnknownsTetra();
    return result;
}

/// \brief Get number of unknowns on a vertex
Uint ParMultiGridCL::NumberOfUnknownsOnVertex()
{
    Uint result=0;
    for (VecDescPCT::const_iterator it(_VecDesc.begin()), end(_VecDesc.end()); it!=end; ++it)
        result +=(*it)->RowIdx->NumUnknownsVertex();
    return result;
}

/// \brief Get number of unknowns on an edge
Uint ParMultiGridCL::NumberOfUnknownsOnEdge()
{
    Uint result=0;
    for (VecDescPCT::const_iterator it(_VecDesc.begin()), end(_VecDesc.end()); it!=end; ++it)
        result += (*it)->RowIdx->NumUnknownsEdge();
    return result;
}

/// \brief Enlarge the Recieve Buffer
void ParMultiGridCL::EnlargeRecieveBuffer()
{
    if (_RecvBuf.size()>0){
        BufferCT tmpBuf(_RecvBuf);
        _RecvBuf.resize(2*_RecvBuf.size());
        std::copy(tmpBuf.begin(), tmpBuf.end(), _RecvBuf.begin());
        Comment("["<<ProcCL::MyRank()<<"]===> Enlarge recieve buffer from "<<tmpBuf.size()<<" to "<<_RecvBuf.size()<<"!"<<std::endl, DebugParallelC);
    }
    else{
        _RecvBuf.resize(1024);
        Comment("["<<ProcCL::MyRank()<<"]===> Create recieve buffer of size "<<_RecvBuf.size()<<"!"<<std::endl, DebugParallelC);
    }
}

/// \brief Gather unknowns on ghost tetras for sending these to master the tetra
int ParMultiGridCL::GatherUnknownsRef (DDD_OBJ obj, void* buf)
/** Get all values on a ghost tetrahedron, that will be deleted, and put these values into the
    given buffer.
*/
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);

    if (tp->GetPrio()!=PrioKilledGhost || !tp->IsMarkedForNoRef())
        return 1;

    TransferUnkT* buffer    = static_cast<TransferUnkT*>(buf);
    Uint          buffer_pos= 0;

    for (VecDescPCT::const_iterator it(_VecDesc.begin()); it!=_VecDesc.end(); ++it)
    {
        const Uint idx=(*it)->RowIdx->GetIdx(); // this should be the index before the refinement algorithm has been started

        if ((*it)->RowIdx->NumUnknownsVertex())
        { // collect unknowns on vertices
            for (TetraCL::const_VertexPIterator sit(tp->GetVertBegin());
                 sit!=tp->GetVertEnd();
                 ++sit,
                 buffer_pos+=(*it)->RowIdx->NumUnknownsVertex())
            {
                if ((*sit)->Unknowns.Exist() && (*sit)->Unknowns.Exist(idx)){    // check for existance
                    for (Uint i=0; i<(*it)->RowIdx->NumUnknownsVertex(); ++i){
                        buffer[buffer_pos+i].val = (*it)->Data[(*sit)->Unknowns(idx)+i];
                        buffer[buffer_pos+i].idx = idx;
                        buffer[buffer_pos+i].mark= true;
                    }
                }
                else{
                    buffer[buffer_pos].mark=false;
                }
            }
        }
        if ((*it)->RowIdx->NumUnknownsEdge())
        { // collect unknowns on edges, if they exists
            for (TetraCL::const_EdgePIterator sit(tp->GetEdgesBegin());
                 sit!=tp->GetEdgesEnd();
                 ++sit,
                buffer_pos+=(*it)->RowIdx->NumUnknownsEdge())
            {
                if ((*sit)->Unknowns.Exist() && (*sit)->Unknowns.Exist(idx))    // check for existance
                    for (Uint i=0; i<(*it)->RowIdx->NumUnknownsVertex(); ++i){
                        buffer[buffer_pos+i].val = (*it)->Data[(*sit)->Unknowns(idx)+i];
                        buffer[buffer_pos+i].idx = idx;
                        buffer[buffer_pos+i].mark= true;
                    }
                else{
                    buffer[buffer_pos].mark=false;
                }
            }
        }
        if ((*it)->RowIdx->NumUnknownsTetra())
        { // collect unknowns on the tetra itselfe
            if (tp->Unknowns.Exist() && tp->Unknowns.Exist(idx))                // check for existance
                for (Uint i=0; i<(*it)->RowIdx->NumUnknownsTetra(); ++i){
                    buffer[buffer_pos+i].val = (*it)->Data[tp->Unknowns(idx)+i];
                    buffer[buffer_pos+i].idx = idx;
                    buffer[buffer_pos+i].mark= true;
                }
            else{
                buffer[buffer_pos].mark=false;
            }
            buffer_pos+=(*it)->RowIdx->NumUnknownsTetra();
        }
    }
    Assert(buffer_pos==NumberOfUnknownsOnTetra(),
           DROPSErrCL("ParMultiGridCL::GatherUnknowns: I haven't check enough places for unknowns!"),
           DebugParallelNumC);
    return 0;
}

/// \brief Scatter unknowns on master tetras, which have been send from killed ghost tetra
int ParMultiGridCL::ScatterUnknownsRef(DDD_OBJ obj, void* buf)
/** This procedure puts all unknowns that can live on subsimplices or the tetrahedron itselfe
    into a recieve buffer, if the unknowns are not known so far
*/
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);

    if (tp->GetPrio()!=PrioMaster)
        return 1;

    TransferUnkT* buffer    = static_cast<TransferUnkT*>(buf);
    Uint          buffer_pos= 0;

    // for all known indices
    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        // index number and unknowns on simplices
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnVert = _VecDesc[index_type]->RowIdx->NumUnknownsVertex(),
                   numUnkOnEdge = _VecDesc[index_type]->RowIdx->NumUnknownsEdge(),
                   numUnkOnTetra= _VecDesc[index_type]->RowIdx->NumUnknownsTetra();

        // Check if there are unknowns on vertices and recieve them
        if (numUnkOnVert)
        {
            // iterate over all vertices of this tetra
            for (TetraCL::const_VertexPIterator sit(tp->GetVertBegin()); sit!=tp->GetVertEnd(); ++sit)
            {
                // put the recieved unknowns into _RecvBuf, if this vertex should store unknowns (i.e. is master)
                // and if the unknown is not known so far.
                if (!(*sit)->Unknowns.Exist(idx)                        // this vert stores no "old" unknown
                     && !(*sit)->Unknowns.UnkRecieved(idx)              // this vert has no unknowns recieved yet
                     && buffer[buffer_pos].mark                         // there have been unknowns on sender side
                     && buffer[buffer_pos].idx==idx                     // and right index
                   )
                {
                    if (_RecvBufPos+numUnkOnVert>_RecvBuf.size())
                        EnlargeRecieveBuffer();
                    (*sit)->Unknowns.Prepare(idx);                      // create UnknownIdxCL
                    (*sit)->Unknowns(idx)= _RecvBufPos;                 // remeber position, where the unknowns has been put

                    for (Uint i=0; i<numUnkOnVert; ++i)                 // store all unknowns
                        _RecvBuf[_RecvBufPos+i]=buffer[buffer_pos+i].val;
                    _RecvBufPos+=numUnkOnVert;

                    (*sit)->Unknowns.SetUnkRecieved(idx);               // this vertex has recieved unknowns
                }
                buffer_pos+=numUnkOnVert;                               // unknowns on vertices has been handled
            }
        }

        // for documentation look at vertices above
        if (numUnkOnEdge)
        {
            for (TetraCL::const_EdgePIterator sit(tp->GetEdgesBegin()); sit!=tp->GetEdgesEnd(); ++sit)
            {
                if ((*sit)->MayStoreUnk()
                     && !(*sit)->Unknowns.Exist(idx)
                     && !(*sit)->Unknowns.UnkRecieved(idx)
                     && buffer[buffer_pos].mark
                     && buffer[buffer_pos].idx==idx
                   )
                {
                    if (_RecvBufPos + numUnkOnEdge>_RecvBuf.size())
                        EnlargeRecieveBuffer();
                    (*sit)->Unknowns.Prepare(idx);
                    (*sit)->Unknowns(idx)= _RecvBufPos;

                    for (Uint i=0; i<numUnkOnEdge; ++i)
                        _RecvBuf[_RecvBufPos+i]=buffer[buffer_pos+i].val;
                    _RecvBufPos+=numUnkOnEdge;

                    (*sit)->Unknowns.SetUnkRecieved(idx);
                }
                buffer_pos+=numUnkOnEdge;
            }
        }

        // Recieve unknowns on tetra itselfe
        if (numUnkOnTetra)
        {
            if (tp->MayStoreUnk()
                    && !tp->Unknowns.Exist(idx)
                    && !tp->Unknowns.UnkRecieved(idx)
                    && buffer[buffer_pos].mark
                    && buffer[buffer_pos].idx==idx
               )
            {
                if (_RecvBufPos + numUnkOnTetra>_RecvBuf.size())
                    EnlargeRecieveBuffer();

                tp->Unknowns.Prepare(idx);
                tp->Unknowns(idx)= _RecvBufPos;

                for (Uint i=0; i<numUnkOnTetra; ++i)
                    _RecvBuf[_RecvBufPos+i] = buffer[buffer_pos+i].val;
                _RecvBufPos+=numUnkOnTetra;

                tp->Unknowns.SetUnkRecieved(idx);
            }
            buffer_pos+=numUnkOnTetra;
        }
    }

    Assert(buffer_pos==NumberOfUnknownsOnTetra(),
           DROPSErrCL("ParMultiGridCL::ScatterUnknowns: I haven't check enough places for unknowns!"),
           DebugParallelNumC);
    return 0;
}

int ParMultiGridCL::GatherUnknownsMigV (DDD_OBJ obj, void* buf)
{
    VertexCL* const sp= ddd_cast<VertexCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;

    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnVert = _VecDesc[index_type]->RowIdx->NumUnknownsVertex();
        if (sp->Unknowns.Exist(idx))
        {
            buffer[bufferpos].mark=true;
            buffer[bufferpos].idx=idx;
            const Uint sysnum= sp->Unknowns(idx);
            if (!sp->Unknowns.UnkRecieved(idx))
                for (Uint i=0; i<numUnkOnVert; ++i)
                    buffer[bufferpos++].val= _VecDesc[index_type]->Data[sysnum+i];
            else //sp->Unknowns.UnkRecieved(idx)
                for (Uint i=0; i<numUnkOnVert; ++i)
                    buffer[bufferpos++].val= _RecvBuf[sysnum+i];
        }
        else
        {
            buffer[bufferpos].mark=false;
            bufferpos+=numUnkOnVert;
        }
    }
    return 0;
}

int ParMultiGridCL::ScatterUnknownsMigV(DDD_OBJ obj, void* buf)
{
    VertexCL* const sp= ddd_cast<VertexCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;

    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnVert = _VecDesc[index_type]->RowIdx->NumUnknownsVertex();
        if (sp->Unknowns.Exist(idx))        // do nothing, unknowns allready known
            bufferpos += numUnkOnVert;
        else if (buffer[bufferpos].mark && buffer[bufferpos].idx==idx)    // new unknowns has been sent to right index
        {
            if (numUnkOnVert+_RecvBufPos>_RecvBuf.size())
                EnlargeRecieveBuffer();

            sp->Unknowns.Prepare(idx);
            sp->Unknowns(idx)=_RecvBufPos;

            for (Uint i=0; i<numUnkOnVert; ++i)
                _RecvBuf[_RecvBufPos+i]=buffer[bufferpos++].val;
            _RecvBufPos+=numUnkOnVert;

            sp->Unknowns.SetUnkRecieved(idx);
        }
        else                                // not known unknowns and no information recieved
            bufferpos+=numUnkOnVert;
    }
    return 0;
}

int ParMultiGridCL::GatherUnknownsMigE (DDD_OBJ obj, void* buf)
{
    EdgeCL* const sp= ddd_cast<EdgeCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;
    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnEdge = _VecDesc[index_type]->RowIdx->NumUnknownsEdge();
        if (!numUnkOnEdge)
            continue;
        if (sp->Unknowns.Exist(idx))
        {
            buffer[bufferpos].mark=true;
            buffer[bufferpos].idx=idx;
            const Uint sysnum= sp->Unknowns(idx);
            if (!sp->Unknowns.UnkRecieved(idx))
                for (Uint i=0; i<numUnkOnEdge; ++i)
                    buffer[bufferpos++].val= _VecDesc[index_type]->Data[sysnum+i];
            else //sp->Unknowns.UnkRecieved(idx)
                for (Uint i=0; i<numUnkOnEdge; ++i)
                    buffer[bufferpos++].val= _RecvBuf[sysnum+i];
        }
        else
        {
            buffer[bufferpos].mark=false;
            buffer[bufferpos].idx=idx;
            bufferpos+=numUnkOnEdge;
        }
    }
    return 0;
}

int ParMultiGridCL::ScatterUnknownsMigE(DDD_OBJ obj, void* buf)
{
    EdgeCL* const sp= ddd_cast<EdgeCL*>(obj);
    TransferUnkT* buffer = static_cast<TransferUnkT*>(buf);
    int bufferpos=0;

    for (size_t index_type=0; index_type<_VecDesc.size(); ++index_type)
    {
        const Uint idx          = _VecDesc[index_type]->RowIdx->GetIdx(),
                   numUnkOnEdge = _VecDesc[index_type]->RowIdx->NumUnknownsEdge();
        if (!numUnkOnEdge)
            continue;
        if (sp->Unknowns.Exist(idx))        // do nothing, unknowns allready known
            bufferpos += numUnkOnEdge;
        else if (buffer[bufferpos].mark && buffer[bufferpos].idx==idx)    // new unknowns has been sent to right index
        {
            if (numUnkOnEdge+_RecvBufPos>_RecvBuf.size())
                EnlargeRecieveBuffer();
            sp->Unknowns.Prepare(idx);
            sp->Unknowns(idx)=_RecvBufPos;
            sp->Unknowns.SetUnkRecieved(idx);

            for (Uint i=0; i<numUnkOnEdge; ++i)
                _RecvBuf[_RecvBufPos+i]=buffer[bufferpos++].val;
            _RecvBufPos+=numUnkOnEdge;
        }
        else                                // not known unknowns and no information recieved
            bufferpos+=numUnkOnEdge;
    }
    return 0;
}


/// \brief Clear all XferNew-Marks
void ParMultiGridCL::DelAllUnkRecv()
{
    for (MultiGridCL::const_VertexIterator sit=_mg->GetAllVertexBegin(); sit!=_mg->GetAllVertexEnd(); ++sit){
        sit->Unknowns.ResetUnkRecieved();
    }
    for (MultiGridCL::const_EdgeIterator sit=_mg->GetAllEdgeBegin(); sit!=_mg->GetAllEdgeEnd(); ++sit){
        sit->Unknowns.ResetUnkRecieved();
    }
    for (MultiGridCL::const_TetraIterator sit=_mg->GetAllTetraBegin(); sit!=_mg->GetAllTetraEnd(); ++sit){
        sit->Unknowns.ResetUnkRecieved();
    }
}

/// \brief Marks all tetras on last level!
void ParMultiGridCL::MarkAll()
{
    for (MultiGridCL::TriangTetraIteratorCL it=_mg->GetTriangTetraBegin(); it!=_mg->GetTriangTetraEnd(); ++it)
        it->SetRegRefMark();
}

/// \brief Refine the Multigrid by calling the procedure ParRefine from the MulitGridCL
void ParMultiGridCL::Refine()
{
    _mg->Refine();
}

/****************************************************************************
* I N I T   A N D   E N D   X F E R                                         *
*****************************************************************************
* Call these functions to start and end xfer-commands                       *
* They are preparing and finishing everything for the transfer.             *
****************************************************************************/
/// \brief prepare everything for the transfers
///
/// Call this everytime before using an transfer command!
void ParMultiGridCL::XferStart(int Level)
{
    Assert(!TransferMode, DROPSErrCL("ParMultiGridCL: XferStart: Allready called XferStart"), DebugParallelC);
    TransferMode=true;          // Now the Transfer-Mode has been started
    ToHandleTetra_.clear();     // There is no tetra that has changed prio from ghost to master within transfer so far

    _mg->ClearTriangCache();    // Remove old iterators, new iterators are comming!
    _mg->PrepareModify();       // the recieved simplices have to be stored, so be prepared for that

    // all procs should have the same number of levels!
    AdjustLevel();

    _level = (Level==-1) ? _mg->GetLastLevel() : Level;

    Comment("- Starting Xfer"<<std::endl, DebugParallelC);
    DDD_XferBegin();
}

/// \brief End the transfer phase by sending simplices, delete unused simplices and eventually rescue unknowns
/** */
void ParMultiGridCL::XferEnd()
{
    Assert(TransferMode && _level!=-1, DROPSErrCL("ParMultiGridCL: XferEnd: Not in Transfer-Mode"), DebugParallelC);
    DDD_XferEnd();

    // All Tetraeders, that are marked for removement, should be removed now!
    // All Subsimplices are marked for removement
    for (int l=0; l<=_level; ++l)
    {
        (*_TetraCont)[l].remove_if( std::mem_fun_ref(&TetraCL::IsMarkedForRemovement) );
        for_each( _mg->GetVerticesBegin(l), _mg->GetVerticesEnd(l), std::mem_fun_ref( &VertexCL::SetRemoveMark ) );
        for_each( _mg->GetEdgesBegin(l), _mg->GetEdgesEnd(l), std::mem_fun_ref( &EdgeCL::SetRemoveMark ) );
        for_each( _mg->GetFacesBegin(l), _mg->GetFacesEnd(l), std::mem_fun_ref( &FaceCL::SetRemoveMark ) );
    }

    DDD_XferBegin();        // for removement/priochange of subsimplices

    Comment("  * Rescue HasGhosts"<<std::endl,DebugParallelC);
    TreatHasGhosts();       // mark subs of HasGhosts as VGhost, rescue subs
    Comment("  * Rescue Ghosts"<<std::endl,DebugParallelC);
    TreatGhosts();          // rescue subs of Ghosts, mark as Ghost
    Comment("  * Rescue Subs"<<std::endl,DebugParallelC);
    // Rescue subs of tetras and remove link to parent, if tetra is ghost
    for (int l=0; l<=_level; ++l){
        for (MultiGridCL::TetraIterator sit((*_TetraCont)[l].begin()), tEnd((*_TetraCont)[l].end()); sit!=tEnd; ++sit){
            RescueSubs(*sit);
            if (sit->IsGhost() && sit->GetParent() ){
                sit->_Parent=0;
            }
        }
    }

    Comment("  * Tell DDD to delete Objects!"<<std::endl,DebugParallelC);
    for (int l=0; l<=_level; ++l)
    {
        for (MultiGridCL::VertexIterator sit((*_VertCont)[l].begin()), tEnd((*_VertCont)[l].end()); sit!=tEnd; ++sit)
        {
            if ( sit->IsMarkedForRemovement() ) sit->XferDelete();
        }
        for (MultiGridCL::EdgeIterator sit((*_EdgeCont)[l].begin()), tEnd((*_EdgeCont)[l].end()); sit!=tEnd; ++sit)
        {
            if ( sit->IsMarkedForRemovement() ) sit->XferDelete();
        }
        for (MultiGridCL::FaceIterator sit((*_FaceCont)[l].begin()), tEnd((*_FaceCont)[l].end()); sit!=tEnd; ++sit)
        {
            if ( sit->IsMarkedForRemovement() ) sit->XferDelete();
        }
    }

    DDD_XferEnd();

    // Adapt midvertex pointers on VGhost-Edges
    Comment("  * Adapting Midvertex on VGhost-Edges"<<std::endl,DebugParallelC);
    AdaptMidVertex();

    // Accumulate Ref-counter on edges
    Comment("  * Accumulate MFR"<<std::endl,DebugParallelC);
    AccumulateMFR();

    // Rescue unknowns on edges, that are deleted and midvertex stays on processor
    if (VecDescRecv()){
        Comment("  * Send unknowns "<<std::endl,DebugParallelC);
        DDD_IFExchange(AllSimplexIFCL<VertexCL>::GetIF(),               // exchange datas over distributed vertices
                       NumberOfUnknownsOnVertex()* sizeof(TransferUnkT),// number of datas to be exchanged
                       GatherUnknownsMigVC,                             // how to gather datas
                       ScatterUnknownsMigVC                             // how to scatter datas
                      );
        DDD_IFExchange(AllSimplexIFCL<EdgeCL>::GetIF(),                 // exchange datas over distributed vertices
                       NumberOfUnknownsOnEdge()* sizeof(TransferUnkT),  // number of datas to be exchanged
                       GatherUnknownsMigEC,                             // how to gather datas
                       ScatterUnknownsMigEC                             // how to scatter datas
                      );
        RescueUnknownsOnEdges();
    }

    // now physically delete the simplices with RemoveMark from memory
    Comment("- Delete all unused Simplices!"<<std::endl,DebugParallelC);
    for (int l=0; l<=_level; ++l)
    {
        (*_VertCont)[l].remove_if( std::mem_fun_ref(&VertexCL::IsMarkedForRemovement) );
        (*_EdgeCont)[l].remove_if( std::mem_fun_ref(&EdgeCL::IsMarkedForRemovement) );
        (*_FaceCont)[l].remove_if( std::mem_fun_ref(&FaceCL::IsMarkedForRemovement) );
    }

    TransferMode = false;       // Now the Transfer-Mode has been ended
    _level=-1;                  // and so _level is also set on no transfer active!
    _mg->FinalizeModify();      // No more elements may be added
    _mg->ClearTriangCache();

    Comment("- Xfer finished"<<std::endl,DebugParallelC);
}


/// \brief This function assures, that every proc has the same number of level. This is important in the refinement algorithm
/** Make sure, the containers are modifyable*/
void ParMultiGridCL::AdjustLevel()
{
    int myLastLevel  =_mg->GetLastLevel();
    int allLastLevel = GlobalMax(myLastLevel);       // this procedure is from parallel.h

    // all procs should have the same number of levels!
    for (; myLastLevel<allLastLevel; ++myLastLevel)
    {
        _mg->AppendLevel();
    }
}

/****************************************************************************
* C H E C K I N G   F O R   S A N I T Y                                          *
*****************************************************************************
*   Functions for checking the sanity of the parallel structures                 *
****************************************************************************/
/// \brief Check for sanity
/** Check if the edges and faces have the same subsimpilces (check with GID).
    Interface-function, calls Gather-, ScatterEdgeSane (EdgeIF)
    and Gather- ScatterFaceSane (FaceIF) and check multigrid for sanity.*/
bool ParMultiGridCL::IsSane(std::ostream& os, int Level)
{
    bool sane= true;

    // Checking all levels
    if (Level==-1)
    {
        Uint maxLevel= GlobalMax( _mg->GetLastLevel());
        if (maxLevel != _mg->GetLastLevel() )
        {
            sane= false;
            os << "Local MultiGrid has too few levels: " << _mg->GetLastLevel()
                    << " instead of " << maxLevel <<std::endl;
        }
        for (Uint lvl=0; lvl<=maxLevel; ++lvl)
        {
            _sane= ParMultiGridCL::IsSane( os, lvl);
            sane= sane && _sane;
        }
    }
    // checking on level
    else
    {
        Comment("Checking level " << Level << std::endl,DebugParallelC);

        sane= _mg->IsSane( os, Level);
        _sane= true; _os= &os;

        Comment("Checking inter-processor dependencies on level " << Level << std::endl, DebugParallelC);

        DDD_IFAOnewayX( _EdgeIF, Level, IF_FORWARD, 2*sizeof(DDD_GID), &GatherEdgeSaneC, &ScatterEdgeSaneC);
        DDD_IFAOnewayX( _FaceIF, Level, IF_FORWARD, 6*sizeof(DDD_GID), &GatherFaceSaneC, &ScatterFaceSaneC);

        sane= sane && _sane;
    }
    return sane;
}

/** Checking if the edges have the same both vertices on all procs that share this edge,
    collect the gids */
int ParMultiGridCL::GatherEdgeSane( DDD_OBJ obj, void* data, DDD_PROC, DDD_ATTR)
{
    DDD_GID* vert= static_cast<DDD_GID*>(data);
    EdgeCL* ep= ddd_cast<EdgeCL*>(obj);
    vert[0]= ep->GetVertex(0)->GetGID();
    vert[1]= ep->GetVertex(1)->GetGID();
    return 0;
}

/** Checking if the edges have the same both vertices on all procs that share this edge,
    recieve the gids */
int ParMultiGridCL::ScatterEdgeSane( DDD_OBJ obj, void* data, DDD_PROC proc, DDD_ATTR)
{
    DDD_GID* vert= static_cast<DDD_GID*>(data);
    EdgeCL* ep= ddd_cast<EdgeCL*>(obj);
    if (vert[0]!=ep->GetVertex(0)->GetGID() || vert[1]!=ep->GetVertex(1)->GetGID())
    {
        *_os << "Vertices of distributed edge differ on proc " << proc
                << " (Verts " << vert[0] << ", " << vert[1] << "). Local edge:\n";
        ep->DebugInfo( *_os);
        _sane= false;
    }
    return 0;
}

/** Checking if the faces have the same three vertices and three edges on all procs that share this edge
    collect information */
int ParMultiGridCL::GatherFaceSane( DDD_OBJ obj, void* data, DDD_PROC, DDD_ATTR)
{
    DDD_GID* vertedge= static_cast<DDD_GID*>(data);
    FaceCL* fp= ddd_cast<FaceCL*>(obj);
    vertedge[0]= fp->GetVertex(0)->GetGID();
    vertedge[1]= fp->GetVertex(1)->GetGID();
    vertedge[2]= fp->GetVertex(2)->GetGID();
    vertedge[3]= fp->GetEdge(0)->GetGID();
    vertedge[4]= fp->GetEdge(1)->GetGID();
    vertedge[5]= fp->GetEdge(2)->GetGID();
    return 0;
}

/** Checking if the faces have the same three vertices and three edges on all procs that share this edge
    recieve information */
int ParMultiGridCL::ScatterFaceSane( DDD_OBJ obj, void* data, DDD_PROC proc, DDD_ATTR)
{
    DDD_GID* ve= static_cast<DDD_GID*>(data);
    FaceCL* fp= ddd_cast<FaceCL*>(obj);
    if (ve[0]!=fp->GetVertex(0)->GetGID() || ve[1]!=fp->GetVertex(1)->GetGID() || ve[2]!=fp->GetVertex(2)->GetGID())
    {
        *_os << "Vertices of distributed face differ on proc " << proc
                << " (Verts " << ve[0] << ", " << ve[1] << ", "<< ve[2] << "). Local face:\n";
        fp->DebugInfo( *_os);
        _sane= false;
    }
    if (ve[3]!=fp->GetEdge(0)->GetGID() ||
           ve[4]!=fp->GetEdge(1)->GetGID() ||
           ve[5]!=fp->GetEdge(2)->GetGID())
    {
        *_os << "Edges of distributed face differ on proc " << proc
                << " (Verts " << ve[0] << ", " << ve[1] << ", "<< ve[2] << "). Local face:\n";
        fp->DebugInfo( *_os);
        _sane= false;
    }
    return 0;
}

/****************************************************************************
* F U N C T I O N S   F O R   I N T E R F A C E S                            *
*****************************************************************************
* The following functions are handlers to call the DDD-Interfaces correct    *
****************************************************************************/
/// \brief accumulate the marks for refinement on Edges on Level (-1==all levels)
/** Interface-function (EdgeIF) calls Gather- and ScatterEdgeMFR
    In the case, a tetra has changed its prio from ghost to master, local edges may
    have inconsistent MFR. This MFR on edges are repaired here, too.*/
void ParMultiGridCL::AccumulateMFR( int Level)
{
    if (Level==-1)
        DDD_IFOneway( _EdgeIF, IF_FORWARD, sizeof(short int), &GatherEdgeMFRC, &ScatterEdgeMFRC);
    else
        DDD_IFAOneway( _EdgeIF, Level, IF_FORWARD, sizeof(short int), &GatherEdgeMFRC, &ScatterEdgeMFRC);

    for (TetraPCT::iterator it(ToHandleTetra_.begin()); it!=ToHandleTetra_.end(); ++it)
        for (TetraCL::const_EdgePIterator epiter((*it)->GetEdgesBegin()); epiter!=(*it)->GetEdgesEnd(); ++epiter)
            if ((*epiter)->IsLocal() && (*epiter)->GetAccMFR()!=(*epiter)->GetMFR())
                (*epiter)->_AccMFR=(*epiter)->_MFR;
    ToHandleTetra_.clear();
}

/** Set AccMFR to MFR and put this information into the message */
int ParMultiGridCL::GatherEdgeMFR( DDD_OBJ obj, void* buf)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);
    short int* buffer= static_cast<short int*>(buf);
    if (!ep->IsMarkedForRemovement())
        *buffer= ep->_AccMFR= ep->_MFR;
    else
        *buffer=0;

    AllComment("- Submit MFR from edge: "<<ep->GetGID()<<" as "<<ep->_AccMFR<<std::endl, DebugParallelHardC);
    return 0;
}

/** increase value of AccMFR */
int ParMultiGridCL::ScatterEdgeMFR( DDD_OBJ obj, void* buf)
{
    EdgeCL* ep=ddd_cast<EdgeCL*>(obj);

    ep->_AccMFR += *static_cast<short int*>(buf);

    AllComment("- Recieve MFR from edge: "<<ep->GetGID()<<", now "<<ep->_AccMFR<<std::endl, DebugParallelHardC);
    return 0;
// Ghost-Kanten enthalten hinterher ebenfalls globalen Wert!
}


/// \brief Communicate the Refinement Marks from all procs
/** Interface-function (TetraIF) calls Gather- and ScatterTetraRestrictMarks */
void ParMultiGridCL::CommunicateRefMarks( Uint Level)
{
    DDD_IFAOneway( _TetraIF, Level, IF_FORWARD, sizeof(bool), &GatherTetraRestrictMarksC, &ScatterTetraRestrictMarksC);
}

/** Put flag, if the Tetra is marked for Refinement into message */
int ParMultiGridCL::GatherTetraRestrictMarks( DDD_OBJ obj, void* buf)
{
    TetraCL* const tp = ddd_cast<TetraCL*>(obj);
    *static_cast<bool*>(buf)= tp->IsMarkedForRegRef();
    AllComment("- Submit RegRefMark from tetra: "<<tp->GetGID()<<" as "<<tp->IsMarkedForRegRef()<<std::endl, DebugParallelHardC);
    return 0;
}

/** Recieve, if the tetra is marked for refinement on an other proc and increas or decrease the MFR on all edges*/
int ParMultiGridCL::ScatterTetraRestrictMarks( DDD_OBJ obj, void* buf)
// Hole Arbeit nach, die in TetraCL::RestrictMarks nicht verrichtet werden konnte.
{
    bool MarkForRegRef= *static_cast<bool*>(buf);
    TetraCL* tp= ddd_cast<TetraCL*>(obj);
    AllComment("- Recieve RegRefMark from tetra: "<<tp->GetGID()<<" as "<<MarkForRegRef<<std::endl,DebugParallelHardC);

    if (tp->IsRegularlyRef() )
    {
        if (!MarkForRegRef)
        {
            tp->SetNoRefMark();
            tp->UnCommitRegRefMark();
        }
    }
    else // tetra is irregularly refined
    {

        Assert( !tp->IsUnrefined(), DROPSErrCL("ParMultiGridCL: ScatterTetraRestrictMarks: Master has Ghost eventhough it is unrefined!"), DebugParallelC);
        if (MarkForRegRef)
        {
            tp->SetRegRefMark();
            tp->CommitRegRefMark();
        }
        else
            tp->SetNoRefMark();
    }
    return 0;
}

/// \brief Treat Ghosts, so they are not deleted
/** Removes the del-mark on all ghost on level k+1 without NoRefMark, so they are not deleted.
    This proecude markes them also as Ghost <p>
    Interface-function (TetraIF) calls ExecGhostRescue */
void ParMultiGridCL::TreatGhosts( int Level)
{
    if (Level==-1)
        DDD_IFExecLocal( _TetraIF, &ExecGhostRescueC);
    else
        DDD_IFAExecLocal( _TetraIF, Level, &ExecGhostRescueC);
}

/** Set Prio for all subsimplices on PrioGhost and removes their Remove-Mark*/
int ParMultiGridCL::ExecGhostRescue( DDD_OBJ obj)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    if (!tp->IsGhost() ) return 1;
    if (!tp->IsMarkedForNoRef())
    {  // rescue all subsimplices and set their Prio to PrioGhost
        for (TetraCL::const_VertexPIterator it= tp->GetVertBegin(), end= tp->GetVertEnd(); it!=end; ++it)
            PrioChange( *it, PrioGhost);
        for (TetraCL::const_EdgePIterator it= tp->GetEdgesBegin(), end= tp->GetEdgesEnd(); it!=end; ++it)
            PrioChange( *it, PrioGhost);
        for (TetraCL::const_FacePIterator it= tp->GetFacesBegin(), end= tp->GetFacesEnd(); it!=end; ++it)
            PrioChange( *it, PrioGhost);
        std::for_each( tp->GetVertBegin(), tp->GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
        std::for_each( tp->GetEdgesBegin(), tp->GetEdgesEnd(), std::mem_fun( &EdgeCL::ClearRemoveMark) );
        std::for_each( tp->GetFacesBegin(), tp->GetFacesEnd(), std::mem_fun( &FaceCL::ClearRemoveMark) );
    }
    return 0;
}

/// \brief Treat Ghost-Vertices, so they are not deleted
/** Vertices has to be treated special for removement, because they can be found in other
    levels then the tetra. <p>
    Interface-function (TetraIF) ExecGhVertRescue */
void ParMultiGridCL::RescueGhostVerts( Uint Level)
{
    // Vertices on  <level> can only owned by ghost tetras on  <Level> to <LastLevel-1>.
    for (Uint lvl= Level, maxLevel= _mg->GetLastLevel(); lvl<maxLevel; ++lvl)
        DDD_IFAExecLocal( _TetraIF, lvl, &ExecGhVertRescueC);
}

/** Rescue vertices of ghost tetras */
int ParMultiGridCL::ExecGhVertRescue( DDD_OBJ obj)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    if (!tp->IsGhost()) return 1;

    if (!tp->IsMarkedForNoRef())
        std::for_each( tp->GetVertBegin(), tp->GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );

    return 0;
}

/// \brief Treat subsimplices that has Ghosts
/** Mark subsimplices of HasGhosts as VGhost and rescue them <p>
    Interface-function (TetraIF) calls ExecHasGhostC*/
void ParMultiGridCL::TreatHasGhosts( int Level)
{
    if (Level==-1)
        DDD_IFExecLocal( _TetraIF, &ExecHasGhostC);
    else
        DDD_IFAExecLocal( _TetraIF, Level, &ExecHasGhostC);
}

/** Set prio for all subsimpilces to PrioVGhost and delete the pointers to children */
int ParMultiGridCL::ExecHasGhost( DDD_OBJ obj)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    if (!tp->HasGhost() ) return 1;

    // This tetra will be removed, so ignore it!
    if (tp->IsGhost() && tp->IsMarkedForNoRef()) return 1;

    std::for_each( tp->GetVertBegin(), tp->GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
    std::for_each( tp->GetEdgesBegin(), tp->GetEdgesEnd(), std::mem_fun( &EdgeCL::ClearRemoveMark) );
    std::for_each( tp->GetFacesBegin(), tp->GetFacesEnd(), std::mem_fun( &FaceCL::ClearRemoveMark) );

    for (TetraCL::const_VertexPIterator it= tp->GetVertBegin(), end= tp->GetVertEnd(); it!=end; ++it)
        PrioChange( *it, PrioVGhost);
    for (TetraCL::const_EdgePIterator it= tp->GetEdgesBegin(), end= tp->GetEdgesEnd(); it!=end; ++it)
        PrioChange( *it, PrioVGhost);
    for (TetraCL::const_FacePIterator it= tp->GetFacesBegin(), end= tp->GetFacesEnd(); it!=end; ++it)
        PrioChange( *it, PrioVGhost);

    if (tp->_Children) { delete tp->_Children; tp->_Children= 0; }
    return 0;
}

/// \brief Adapt midvertex pointers on VGhost-Edges
/** Interface-function (EdgeIF) calls ExecAdaptVGhostMidVertexC*/
void ParMultiGridCL::AdaptMidVertex( int Level)
{
    if (Level==-1)
        DDD_IFExecLocal( _EdgeIF, &ExecAdaptVGhostMidVertexC);
    else
        DDD_IFAExecLocal( _EdgeIF, Level, &ExecAdaptVGhostMidVertexC);
}

/** Delete MidVertex of PrioVGhost Kanten*/
int ParMultiGridCL::ExecAdaptVGhostMidVertex( DDD_OBJ obj)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);

    if (ep->GetPrio()!=PrioVGhost) return 1;

    if (ep->GetMidVertex() && ep->GetMidVertex()->IsMarkedForRemovement())
    {
        ep->RemoveMidVertex();
    }
    return 0;
}

/// \brief Set prios of all subsimplices right
/** This uses TreatHasGhosts(), TreatGhosts() and RescueSubs() */
void ParMultiGridCL::AdaptPrioOnSubs()
{
    DDD_XferBegin();  // for removement/priochange of subsimplices
    TreatHasGhosts(); // mark subs of HasGhosts as VGhost, rescue subs
    TreatGhosts();    // rescue subs of Ghosts, mark as Ghost
    for (Uint level=0; level<=_mg->GetLastLevel(); ++level)
    {
        for (MultiGridCL::TetraIterator sit((*_TetraCont)[level].begin()), tEnd((*_TetraCont)[level].end()); sit!=tEnd; ++sit)
            if (!sit->IsGhost() || !sit->IsMarkedForNoRef())
                RescueSubs(*sit);
        //for_each( _mg->GetTetrasBegin(level), _mg->GetTetrasEnd(level), RescueSubs);
    }
    DDD_XferEnd();
}

/// \brief Rescue all subsimplices and set prios to PrioMaster
void ParMultiGridCL::RescueSubs( TetraCL& t)
{
    if (t.IsMarkedForRemovement() || t.HasGhost() || t.IsGhost() ) return;
    // keep my sub simplices!
    std::for_each( t.GetVertBegin(), t.GetVertEnd(), std::mem_fun( &VertexCL::ClearRemoveMark) );
    std::for_each( t.GetEdgesBegin(), t.GetEdgesEnd(), std::mem_fun( &EdgeCL::ClearRemoveMark) );
    std::for_each( t.GetFacesBegin(), t.GetFacesEnd(), std::mem_fun( &FaceCL::ClearRemoveMark) );

        // mark my subs with prio Master
    for (TetraCL::const_VertexPIterator it= t.GetVertBegin(), end= t.GetVertEnd(); it!=end; ++it)
        PrioChange( *it, PrioMaster);
    for (TetraCL::const_EdgePIterator it= t.GetEdgesBegin(), end= t.GetEdgesEnd(); it!=end; ++it)
        PrioChange( *it, PrioMaster);
    for (TetraCL::const_FacePIterator it= t.GetFacesBegin(), end= t.GetFacesEnd(); it!=end; ++it)
        PrioChange( *it, PrioMaster);
}

void ParMultiGridCL::MarkSimplicesForUnknowns(int lvl)
{
    DDD_XferBegin();
    for (MultiGridCL::TriangTetraIteratorCL  tit(_mg->GetTriangTetraBegin(lvl)), tend(_mg->GetTriangTetraEnd(lvl));
         tit!=tend; ++tit)
    {
        // master tetras in last triang level are able to store unknowns the rest isn't
        for (TetraCL::const_VertexPIterator it(tit->GetVertBegin()), end(tit->GetVertEnd()); it!=end; ++it)
        {
            Assert((*it)->GetPrio()!=PrioVGhost, DROPSErrCL("ParMultiGridCL::MarkSimplicesForUnknowns: Marking PrioVGhost as PrioHasUnk"), DebugParallelNumC);
            PrioChange( *it, PrioHasUnk);
        }
        for (TetraCL::const_EdgePIterator it(tit->GetEdgesBegin()), end(tit->GetEdgesEnd()); it!=end; ++it)
        {
            Assert((*it)->GetPrio()!=PrioVGhost, DROPSErrCL("ParMultiGridCL::MarkSimplicesForUnknowns: Marking PrioVGhost as PrioHasUnk"), DebugParallelNumC);
            PrioChange( *it, PrioHasUnk);
        }
        for (TetraCL::const_FacePIterator it(tit->GetFacesBegin()), end(tit->GetFacesEnd()); it!=end; ++it)
        {
            Assert((*it)->GetPrio()!=PrioVGhost, DROPSErrCL("ParMultiGridCL::MarkSimplicesForUnknowns: Marking PrioVGhost as PrioHasUnk"), DebugParallelNumC);
            PrioChange( *it, PrioHasUnk);
        }
    }
    DDD_XferEnd();
}

/// \brief Destroy unknowns on non-master vertices, edges and tetras if there are information about them
void ParMultiGridCL::DeleteUnksOnGhosts(int Level)
/** This procedure iterates over all simplices, that stores unknowns and delete them.
    \param Level level of the simplices (default all levels)
*/
{
    if (_UnkOnSimplex[0])
        for (MultiGridCL::VertexIterator sit(_mg->GetAllVertexBegin(Level)), end(_mg->GetAllVertexEnd(Level)); sit!=end; ++sit)
            if (!sit->IsMaster() && sit->Unknowns.Exist())
                sit->Unknowns.Destroy();
    if (_UnkOnSimplex[1])
        for (MultiGridCL::EdgeIterator sit(_mg->GetAllEdgeBegin(Level)), end(_mg->GetAllEdgeEnd(Level)); sit!=end; ++sit)
            if (!sit->IsMaster() && sit->Unknowns.Exist())
                sit->Unknowns.Destroy();
    if (_UnkOnSimplex[2])
        for (MultiGridCL::TetraIterator sit(_mg->GetAllTetraBegin(Level)), end(_mg->GetAllTetraEnd(Level)); sit!=end; ++sit)
            if (!sit->IsMaster() && sit->Unknowns.Exist())
                sit->Unknowns.Destroy();


/*    if (_VecDescRecv){
        if (Level==-1)
        {
            DDD_IFExecLocal(NotMasterIF_, &ExecDestroyUnksVC);
            DDD_IFExecLocal(NotMasterIF_, &ExecDestroyUnksEC);
        }
        else
        {
            DDD_IFAExecLocal(NotMasterIF_, Level, &ExecDestroyUnksVC);
            DDD_IFAExecLocal(NotMasterIF_, Level, &ExecDestroyUnksEC);
        }
    }*/
}

template<class SimplexT>
  int ParMultiGridCL::DestroyUnksOnSimplex(DDD_OBJ obj)
{
    SimplexT *sp= ddd_cast<SimplexT*>(obj);

    if (sp->Unknowns.Exist())
    {
        sp->Unknowns.Destroy();
        return 0;
    }
    else
        return 1;
}

/****************************************************************************
* I D E N T I F Y  -  F U N C T I O N S                                          *
*****************************************************************************
*   The folowing functions indentify Vertices, Edges or Faces. If procs      *
*   creates the same subsimplices, this have to be done.                         *
****************************************************************************/
/// \brief Identify a vertex by parent edge
void ParMultiGridCL::IdentifyVertex( const EdgeCL* Parent)
{
    DDD_HDR parhdr= const_cast<DDD_HDR>( Parent->GetHdr()),
    hdr= const_cast<DDD_HDR>( Parent->GetMidVertex()->GetHdr());
    for( int* proclist= DDD_InfoProcList( parhdr)+2; *proclist!=-1; proclist+= 2)
        if (proclist[1]!=PrioVGhost)
        {
            DDD_IdentifyObject( hdr, *proclist, parhdr);
        }
}

/// \brief Identify an edge by parent edge
void ParMultiGridCL::IdentifyEdge( EdgeCL* Me, const EdgeCL* Parent, Uint nr)
{
    DDD_HDR parhdr= const_cast<DDD_HDR>( Parent->GetHdr());
    for( int* proclist= DDD_InfoProcList( parhdr)+2; *proclist!=-1; proclist+= 2)
        if (proclist[1]!=PrioVGhost)
        {
            DDD_IdentifyObject( Me->GetHdr(), *proclist, parhdr);
            DDD_IdentifyNumber( Me->GetHdr(), *proclist, nr);
        }
}

/// \brief Identify an edge by parent face and two vertices
void ParMultiGridCL::IdentifyEdge( EdgeCL* Me, const FaceCL* Parent, const VertexCL* vp0, const VertexCL* vp1)
{
    DDD_HDR parhdr= const_cast<DDD_HDR>( Parent->GetHdr()),
    hdr0= const_cast<DDD_HDR>( vp0->GetHdr()),
    hdr1= const_cast<DDD_HDR>( vp1->GetHdr());
    for( int* proclist= DDD_InfoProcList( parhdr)+2; *proclist!=-1; proclist+= 2)
        if (proclist[1]!=PrioVGhost)
        {
            DDD_IdentifyObject( Me->GetHdr(), *proclist, hdr0);
            DDD_IdentifyObject( Me->GetHdr(), *proclist, hdr1);
        }
}

/// \brief Identify a face with parent face and a number
void ParMultiGridCL::IdentifyFace( FaceCL* Me, const FaceCL* Parent, Uint nr)
{
    DDD_HDR parhdr= const_cast<DDD_HDR>( Parent->GetHdr());
    for( int* proclist= DDD_InfoProcList( parhdr)+2; *proclist!=-1; proclist+= 2)
        if (proclist[1]!=PrioVGhost)
        {
            DDD_IdentifyObject( Me->GetHdr(), *proclist, parhdr);
            DDD_IdentifyNumber( Me->GetHdr(), *proclist, nr);
        }
}

/****************************************************************************
* D E C L A R E A L L                                                            *
*****************************************************************************
*   This procedure declares the DDD_TYPEs: Vertex, Edge, Face, Tetra,        *
*   AddedScal, AddedVec, BndPtr, ChildPtr                                        *
****************************************************************************/
/// \brief Declare all DDD-Types
void ParMultiGridCL::DeclareAll()
{
    VertexCL::Declare();
    EdgeCL::Declare();
    FaceCL::Declare();
    TetraCL::Declare();
    AddedScalCL::Declare();
    AddedVecCL::Declare();
    DeclareBndPtT();
    DeclareChildPtrT();
    Comment("- All Types are declared" << std::endl, DebugParallelC);
}

/****************************************************************************
* D E F I N E A L L                                                             *
*****************************************************************************
*   This procedure defines the DDD_TYPEs: Vertex, Edge, Face, Tetra,         *
*   AddedScal, AddedVec, BndPtr, ChildPtr                                        *
****************************************************************************/
/// \brief Define all DDD-Types
void ParMultiGridCL::DefineAll()
{
    VertexCL::Define();
    EdgeCL::Define();
    FaceCL::Define();
    TetraCL::Define();
    AddedScalCL::Define();
    AddedVecCL::Define();
    DefineBndPtT();
    DefineChildPtrT();
    Comment("- All Types are defined" << std::endl, DebugParallelC);
}

/****************************************************************************
* I N I T  I F                                                                      *
*****************************************************************************
* This procedure definies the Interfaces for: Edges, Faces, Tetras          *
* The interfaces for numerical accumulations are also defined here. So the  *
* user does not have to worry about that!                                       *
****************************************************************************/
/// \brief Initialize all Interfaces
void ParMultiGridCL::InitIF()
{
    DDD_TYPE  O[4];
    O[0]= VertexCL::GetType();
    O[1]= EdgeCL::GetType();
    O[2]= TetraCL::GetType();
    O[3]= FaceCL::GetType();

    DDD_PRIO  A[5], B[4];   // arrays of priorities

    A[0]= B[0] = PrioHasUnk;
    A[1]= B[1] = PrioMaster;
    A[2]= B[2] = PrioVGhost;
    A[3]= B[3] = PrioGhost;
    A[4]=        PrioKilledGhost;

    // interface of edges
    Assert(!_EdgeIF, DROPSErrCL("ParMultiGridCL: InitIF: EdgeIF allready declared"), DebugParallelC);
    _EdgeIF = DDD_IFDefine(1, O+1, 4, A, 4, B);               // Master, VGhost, Ghost -> Master, VGhost, Ghost
    DDD_IFSetName( _EdgeIF, (char*)"Edge-IF for AccumulateMFR");

    // interface of faces
    Assert(!_FaceIF, DROPSErrCL("ParMultiGridCL: InitIF: FaceIF allready declared"), DebugParallelC);
    _FaceIF = DDD_IFDefine(1, O+3, 2, A, 2, B);               // Master Face -> Master Face
    DDD_IFSetName( _FaceIF, (char*)"Face-IF for Sanity Check");

    // interface of tetras
    Assert(!_TetraIF, DROPSErrCL("ParMultiGridCL: InitIF: TetraIF allready declared"), DebugParallelC);
    _TetraIF = DDD_IFDefine(1, O+2, 1, A+3, 2, B);    // Ghost Tetra -> Master Tetra
    DDD_IFSetName( _TetraIF, (char*)"Tetra-IF for RestrictMarks");

    // interface of not master vertices, edges and tetras
    Assert(!NotMasterIF_, DROPSErrCL("ParMultiGridCL: InitIF: NotMasterV_IF allready declared"), DebugParallelC);
    NotMasterIF_= DDD_IFDefine(3, O, 2, A+2, 2, B+2);   // PrioVGhost, PrioGhost -> PrioVGhost, PrioGhost
    DDD_IFSetName( NotMasterIF_, (char*)"non master Vertex, Edge, Tetrahedron-IF");

    // Ghost tetras to master-tetras
    _GhToMaTetra_IF= DDD_IFDefine(1, O+2, 1, A+4, 2, B);    // PrioKilledGhost -> PrioMaster
    DDD_IFSetName( _GhToMaTetra_IF, (char*)"Killed Ghost to Master Interface");

    // Init the Interfaces for numerical accumulations!
    InterfaceCL<VertexCL>::InitIF();
    InterfaceCL<EdgeCL>::InitIF();
    InterfaceCL<TetraCL>::InitIF();
    AllSimplexIFCL<VertexCL>::InitIF();
    AllSimplexIFCL<EdgeCL>::InitIF();
}

/****************************************************************************
* < S i m p l e x > X f e r                                                 *
*****************************************************************************
* These procedures transfer a simplex from the calling proc to another      *
* proc. XferStart() and XferEnd() have to call before and after respectively*
****************************************************************************/
/// \brief Send a Vertex
void ParMultiGridCL::VXfer(VertexCL &v, DDD_PROC dest, DDD_PRIO prio, bool del)
{
    // Set the removement-mark before calling the transfer command, because this mark is checked, whether additional data is send too!
    if (del) v.SetRemoveMark();
    const DDD_HDR hdr= const_cast<DDD_HDR>(&v._dddH);
    DDD_XferCopyObj( hdr, dest, prio);
}

/// \brief Send a Edge
void ParMultiGridCL::EXfer(EdgeCL &e, DDD_PROC dest, DDD_PRIO prio, bool del)
{
    if (del)
        e.SetRemoveMark();

    const DDD_HDR hdr= const_cast<DDD_HDR>(&e._dddH);
    DDD_XferCopyObj( hdr, dest, prio);

}

/// \brief Send a Face
void ParMultiGridCL::FXfer(FaceCL &f, DDD_PROC dest, DDD_PRIO prio, bool del)
{
    if (del)
        f.SetRemoveMark();

    const DDD_HDR hdr= const_cast<DDD_HDR>(&f._dddH);
    DDD_XferCopyObj(hdr, dest, prio);
    // del richtig berücksichtigen!
}

/// \brief Send a Tetra
void ParMultiGridCL::TXfer(TetraCL &t, DDD_PROC dest, DDD_PRIO prio, bool del)
{
    const DDD_HDR hdr= const_cast<DDD_HDR>(&t._dddH);

    DDD_XferCopyObj( hdr, dest, prio);

    if (del)  // not Ma->Gh-Copy
    {
        t.UnlinkFromFaces();
        t.XferDelete();
    }
    if (t.IsRegularlyRef() && t.IsMaster() && prio==PrioMaster)
        t.UnCommitRegRefMark();
}


/****************************************************************************
* H A N D L E R S   F O R   D D D                                           *
*****************************************************************************
* The following handlers tell the DDD-System how to tread the DROPS classes *
* if they are touched by Xfer-Commands or by Interface-Commands             *
****************************************************************************/

//  Deleting an object by DDD
//----------------------------
/// \brief Handle delete of an object
/** Set the remove mark and desctuct the DDD-Hdr. The real delete of the object is done within
    the code.*/
template<class SimplexT> void ParMultiGridCL::HandlerDelete( DDD_OBJ obj)
{
    SimplexT* const sp= ddd_cast<SimplexT*>(obj);
    sp->SetRemoveMark();
    //  Ich denke, dass man hier den Hdr_Destructor nicht braucht, da er automatisch von DDD aufgerufen wird
    // ---> Doch, den braucht man ganz ganz unbedingt! (siehe Hdr_Destructor bei DDD!)
    DDD_HdrDestructor( &sp->_dddH);
    AllComment("  * Simplex with GID "<<sp->GetGID()<<" deleted on proc "<<ProcCL::MyRank()<<std::endl,DebugParallelHardC);
}

//  Construct an object by DDD
//-----------------------------
/// \brief Construct a vertex and return it as a DDD-Object
DDD_OBJ ParMultiGridCL::HandlerVConstructor( size_t, DDD_PRIO, DDD_ATTR level){
    (*_VertCont)[level].push_back( VertexCL());
    return ddd_cast(&(*_VertCont)[level].back());
}

/// \brief Construct a edge and return it as a DDD-Object
DDD_OBJ ParMultiGridCL::HandlerEConstructor( size_t, DDD_PRIO, DDD_ATTR level){
    (*_EdgeCont)[level].push_back( EdgeCL());
    return ddd_cast(&(*_EdgeCont)[level].back());
}

/// \brief Construct a face and return it as a DDD-Object
DDD_OBJ ParMultiGridCL::HandlerFConstructor( size_t, DDD_PRIO, DDD_ATTR level){
    (*_FaceCont)[level].push_back( FaceCL());
    return ddd_cast(&(*_FaceCont)[level].back());
}

/// \brief Construct a tera and return it as a DDD-Object
DDD_OBJ ParMultiGridCL::HandlerTConstructor( size_t, DDD_PRIO, DDD_ATTR level){
    (*_TetraCont)[level].push_back( TetraCL());
    return ddd_cast(&(*_TetraCont)[level].back());
}

//  Transfer an object by DDD
//----------------------------
/// \brief transfer a vertex
/** Check if numerical data have to be transfered too. If this case happens, tell DDD that additional
    data will be transfered. Than DDD calls the Gather and Scatter functions <p>
    The Recylce-Bin is also destroyed and boundary information are send too*/
void ParMultiGridCL::HandlerVXfer( DDD_OBJ obj, __UNUSED__ DDD_PROC proc, __UNUSED__ DDD_PRIO prio)
{
    VertexCL* const vp= ddd_cast<VertexCL*>(obj);
    int numSendScalUnk= 0,
        numSendVecUnk = 0;

    vp->DestroyRecycleBin();

    // if there are bondary-information, then send them too!
    if (vp->IsOnBoundary() )
        DDD_XferAddData( vp->_BndVerts->size(), _BndPtT);

    // if there this ParMultiGridCL knowns about unknowns, ther are unknwons on this vertex and this vertex is marked for removement, count the unknowns and give this Information to DDD!
    if (VecDescRecv() && vp->Unknowns.Exist())
    {
        for (size_t i=0; i<_VecDesc.size(); ++i)
        {
            if ((_VecDesc[i]) && (_VecDesc[i])->RowIdx &&  vp->Unknowns.Exist((_VecDesc[i])->RowIdx->GetIdx()) )
            {
                Assert((_VecDesc[i])->RowIdx->NumUnknownsVertex() ==1 || (_VecDesc[i])->RowIdx->NumUnknownsVertex()==3,
                        DROPSErrCL("ParMultiGridCL::HandlerVXfer: Can only send scalar or 3d-vectors"), DebugParallelC);
                if ( (_VecDesc[i])->RowIdx->NumUnknownsVertex() ==1 )
                    ++numSendScalUnk;
                if( (_VecDesc[i])->RowIdx->NumUnknownsVertex() ==3 )
                    ++numSendVecUnk;
            }
        }
        DDD_XferAddData( numSendScalUnk, AddedScalCL::_dddT );
        DDD_XferAddData( numSendVecUnk,  AddedVecCL::_dddT  );
    }

    AllComment("  * Vertex with GID " << vp->GetGID() << " from " << ProcCL::MyRank() << " to " << proc << " "
            << "with " << ( !(vp->IsOnBoundary()) ? 0 : vp->_BndVerts->size()) << " boundary points and "
            << "with " << numSendScalUnk << " scalar unknowns and "
            << "with " << numSendVecUnk  << " vectoriel unknowns!" << std::endl, DebugParallelHardC);
}

/// \brief transfer an edge
/** Check if numerical data have to be transfered too. If this case happens, tell DDD that additional
    data will be transfered. Than DDD calls the Gather and Scatter functions <p>*/
void ParMultiGridCL::HandlerEXfer(DDD_OBJ obj, __UNUSED__ DDD_PROC proc, DDD_PRIO)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);

    int numSendScalUnk= 0,
        numSendVecUnk = 0;

    // if ParMultiGridCL knowns about unknowns and there are unknowns on this
    // edge count the unknowns and give this information to DDD
    if (VecDescRecv() && ep->Unknowns.Exist())
    {
        for (size_t i=0; i<_VecDesc.size(); ++i)
        {
            if ((_VecDesc[i]) && (_VecDesc[i])->RowIdx &&  ep->Unknowns.Exist((_VecDesc[i])->RowIdx->GetIdx()) )
            {
                Assert((_VecDesc[i])->RowIdx->NumUnknownsEdge()==1 || (_VecDesc[i])->RowIdx->NumUnknownsEdge()==3,
                        DROPSErrCL("ParMultiGridCL::HandlerVXfer: Can only send scalar or 3d-vectors"), DebugParallelC);

                if ( (_VecDesc[i])->RowIdx->NumUnknownsEdge() ==1 )
                    ++numSendScalUnk;
                if ( (_VecDesc[i])->RowIdx->NumUnknownsEdge() ==3 )
                            ++numSendVecUnk;
            }
        }
        DDD_XferAddData( numSendScalUnk, AddedScalCL::_dddT );
        DDD_XferAddData( numSendVecUnk,  AddedVecCL::_dddT  );
    }

    AllComment("  * Edge with GID " << ep->GetGID() << " from " << ProcCL::MyRank() << " to " << proc << " "
            << "with " << numSendScalUnk << " scalar unknowns and "
            << "with " << numSendVecUnk  << " vectoriel unknowns!" << std::endl, DebugParallelHardC);
}

/// \brief transfer a face
/** Nothing is done. <p>*/
void ParMultiGridCL::HandlerFXfer(__UNUSED__ DDD_OBJ obj, __UNUSED__ DDD_PROC proc, DDD_PRIO)
{
    AllComment("  * Face with GID " << ddd_cast<FaceCL*>(obj)->GetGID() << " from " << ProcCL::MyRank() << " to " << proc << std::endl, DebugParallelHardC);
}

/// \brief transfer a tetra
/** Transfer a tetraeder with all subsimplices. Also pointer to childern are transfered and
    if neccessary numerical data.*/
void ParMultiGridCL::HandlerTXfer(DDD_OBJ obj, DDD_PROC proc, DDD_PRIO prio)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    int numChilds=0;
    if (tp->_Children && !tp->HasGhost() )
    {
        DDD_XferAddData( tp->GetRefData().ChildNum, _ChildPtrT);
        numChilds = tp->GetRefData().ChildNum;
    }

    int numSendScalUnk= 0,
        numSendVecUnk = 0;

    // if there this ParMultiGridCL knowns about unknowns, ther are unknwons on this tetra, count the unknowns and give this Information to DDD!
    if (VecDescRecv() && tp->Unknowns.Exist())
    {
        for (size_t i=0; i<_VecDesc.size(); ++i)
        {
            if ((_VecDesc[i]) && (_VecDesc[i])->RowIdx &&  tp->Unknowns.Exist((_VecDesc[i])->RowIdx->GetIdx()) )
            {
                Assert((_VecDesc[i])->RowIdx->NumUnknownsTetra() ==1 || (_VecDesc[i])->RowIdx->NumUnknownsTetra()==3,
                        DROPSErrCL("ParMultiGridCL::HandlerTXfer: Can only send scalar or 3d-vectors"), DebugParallelC);

                if ( (_VecDesc[i])->RowIdx->NumUnknownsTetra()==1 )
                    ++numSendScalUnk;
                if ( (_VecDesc[i])->RowIdx->NumUnknownsTetra()==3 )
                    ++numSendVecUnk;
            }
        }
        DDD_XferAddData( numSendScalUnk, AddedScalCL::_dddT );
        DDD_XferAddData( numSendVecUnk,  AddedVecCL::_dddT  );
    }

    // all subsimplices and of the tetra are transferred too.
    for (Uint i=0; i<NumVertsC; ++i)
        DDD_XferCopyObj( tp->_Vertices[i]->GetHdr(), proc, prio);
    for (Uint i=0; i<NumEdgesC; ++i)
        DDD_XferCopyObj( tp->_Edges[i]->GetHdr(), proc, prio);
    for (Uint i=0; i<NumFacesC; ++i)
        DDD_XferCopyObj( tp->_Faces[i]->GetHdr(), proc, prio);

    AllComment("  * Tetra with GID " << tp->GetGID() << " from " << ProcCL::MyRank() << " to " << proc
            << " with " << numSendScalUnk << " scalar unknowns,"
            << " with " << numSendVecUnk << " vectoriel unknowns"
            << " and " << numChilds << " Children" << std::endl, DebugParallelHardC);
}

//  Gather and Scatter-Handlers
//-----------------------------

/// \brief Send additional data with the simplices
/** These procedures put additional data to a message or recieve this data.
    This data might be geometrical or numerical, like boundary-information or children-information,
    or the numerical values onto the simplex. */
void ParMultiGridCL::HandlerVGather( DDD_OBJ obj, int cnt, DDD_TYPE type, void* buf)
{
    VertexCL* const vp= ddd_cast<VertexCL*>(obj);
    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType() || type==_BndPtT,
           DROPSErrCL("ParMultiGridCL: HandlerVGather: Cannot handle this type!"), DebugParallelC);
    // if this vert is on a boundary, add the boundary-points to the message
    if (type == _BndPtT)
    {
        // the buffer is a storage of boundary-pointers
        BndPointCL* buffer= static_cast<BndPointCL*>(buf);
        // put all boundary vertices into the buffer
        for( VertexCL::const_BndVertIt it= vp->GetBndVertBegin(), end= vp->GetBndVertEnd(); it!=end; ++it, ++buffer)
        {
            //std::cerr << it->GetCoord2D() << std::endl;
            *buffer= *it;
        }
    }
    // if there are Unknowns on this vert, add these to the message
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type == AddedVecCL::GetType()) )
        SendUnknowns(vp,type,buf,cnt);
}

void ParMultiGridCL::HandlerVScatter( DDD_OBJ obj, int cnt, DDD_TYPE type, void* buf, int /*newness*/)
{
    VertexCL* const vp= ddd_cast<VertexCL*>(obj);
    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType() || type==_BndPtT ,
           DROPSErrCL("ParMultiGridCL: HandlerVScatter: Cannot handle this type!"), DebugParallelC);
    // if boundary information are recieved
    if ( type == _BndPtT )
    {
        const BndPointCL* const buffer= static_cast<BndPointCL*>(buf);
//         vp->_BndVerts= new std::vector<BndPointCL>;         // allocate memory for the boundary-points
        for( int i=0; i<cnt; ++i)
            if (!vp->HasBnd(buffer[i]))
                vp->AddBnd(buffer[i]);
//             vp->_BndVerts->push_back( buffer[i]);           // store the recieved boundary-points
    }

    // if numerical data are recieved
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type==AddedVecCL::GetType()) )
        RecvUnknowns(vp,type,buf,cnt);
}

void ParMultiGridCL::HandlerEGather( DDD_OBJ obj, int cnt, DDD_TYPE type, void* buf)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);

    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerEGather: Cannot handle this type!"), DebugParallelC);
    if ( VecDescRecv() && (type == AddedScalCL::GetType() || type == AddedVecCL::GetType()) )
        SendUnknowns(ep,type,buf,cnt);
}

void ParMultiGridCL::HandlerEScatter( DDD_OBJ obj, int cnt, DDD_TYPE type, void* buf, int /*newness*/)
{
    EdgeCL* const ep= ddd_cast<EdgeCL*>(obj);
    Assert(type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerEScatter: Cannot handle this type!"), DebugParallelC);

    if ( VecDescRecv() && (type == AddedScalCL::GetType() || type==AddedVecCL::GetType()) )
        RecvUnknowns(ep,type,buf,cnt);
}

/// \brief Collect additional data for a tetra transfer
/** The additional data may be pointer to children or numerical data within a migration.*/
void ParMultiGridCL::HandlerTGather( DDD_OBJ obj, int cnt, DDD_TYPE type, void* buf)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    Assert(type==_ChildPtrT || type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerTGather: Cannot handle this type!"), DebugParallelC);
    if (type==_ChildPtrT)           // if this Tetra should send his children
    {
        TetraCL** buffer= static_cast<TetraCL**>(buf);
        for (TetraCL::ChildPIterator it(tp->GetChildBegin()), end(tp->GetChildEnd()); it!=end; ++it, ++buffer)
        {
            *buffer= *it;
        }
        // For Problem prob2_sun_2_procs.param (Test-Case for pointer-length in DDD!)
//      if (tp->GetGID()==1728){
//          std::cerr << "["<<ProcCL::MyRank()<<"]  sizeof(*buffer)="<<sizeof(*buffer)<<", sizeof(*ChildPIterator)="<<sizeof(*(tp->GetChildBegin()))<<"  Tetra: "<<tp->GetGID()<<", werde folgende Kinder versenden:\n   ";
//          for (int i=0; i<cnt; ++i)
//              std::cerr << (tmp[i])->GetGID() <<",  ";
//          std::cerr  << std::endl;
//      }
    }
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type == AddedVecCL::GetType()) )
        SendUnknowns(tp,type,buf,cnt);
}

/// \brief Recieve additional data for a tetra transfer
void ParMultiGridCL::HandlerTScatter( DDD_OBJ obj, int cnt, DDD_TYPE type, void* buf, int newness)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    Assert(type==_ChildPtrT || type==AddedScalCL::GetType() || type==AddedVecCL::GetType(),
           DROPSErrCL("ParMultiGridCL: HandlerTScatter: Cannot handle this type!"), DebugParallelC);
    if (newness!=XFER_NEW && tp->IsGhost() )
    {
        // This case shouldn't happen, because this case is catched by the migrate function in LoadBalCL
        std::cerr << ">>>>["<<ProcCL::MyRank()<<"] HandlerScatterTetra: ("<<tp->GetGID()<<") Master replaced by Ghost, continuing anyway, workaround enabled! Seems to be an illegal xfer!!!!"<<std::endl;
        tp->GetHdr()->prio= PrioMaster;
    }

    if (type==_ChildPtrT)                           // hiphiphurra, I am mother!
    {
        TetraCL** const buffer= static_cast<TetraCL**>(buf);

        // For Problem Problems/prob2_sun_2_procs.param! (Test-Case for pointer-length in DDD!)
//      if (ProcCL::MyRank()==0 && tp->GetGID()==1728){
//          std::cerr << "["<<ProcCL::MyRank()<<"]  sizeof(*buffer)="<<sizeof(*buffer)<<", sizeof(buffer)="<<sizeof(buffer)<<": Tetra "<<tp->GetGID()<<" Soll folgende Kinder empfangen:\n   ";
//          for (int i=0; i<cnt; ++i)
//              std::cerr << (tmp[i])->GetGID() << ",  ";
//          std::cerr  << std::endl;
//      }

        // Create new children-container if necessary
        if (!tp->_Children)
            tp->_Children= new SArrayCL<TetraCL*, MaxChildrenC>;
        // put recieved children into children-container!
        for( int i=0; i<cnt; ++i)
            (*(tp->_Children))[i]= buffer[i];
    }
    else if ( VecDescRecv() && (type == AddedScalCL::GetType() || type==AddedVecCL::GetType()) )
        RecvUnknowns(tp,type,buf,cnt);
}

/// \brief Make Edge MFR consistent
/** If a tetra is regular refined and the priority is master and the tetraeder is made on this proc within the actual
    transfer mode, then this tetra increas MFR and set AccMFR to MFR.
    XferEnd() will accumulate the MFRs.*/
void ParMultiGridCL::HandlerTObjMkCons( DDD_OBJ obj, int newness)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);

    if (tp->IsRegularlyRef() && tp->IsMaster() && newness==XFER_NEW)
    {
        tp->CommitRegRefMark();
        // nun wird auf allen Edges _AccMFR:=_MFR gesetzt, um Unkonsistenzen bei vorher verteilt
        // und nun nur noch lokal gespeicherten Edges zu vermeiden. Ein abschliessenden
        // AccumulateMFR in XferEnd() setzt auf den verteilt gespeicherten Edges dann die
        // richtigen _AccMFR-Werte.
        for (TetraCL::const_EdgePIterator it= tp->GetEdgesBegin(), end= tp->GetEdgesEnd(); it!=end; ++it)
            (*it)->_AccMFR= (*it)->_MFR;
    }
}



/// \brief Set Prio of an tetra
/** if priority is set to ghost, the pointer to parent is set to null
    If Prio changes from Ghost to Master, the MFR-Counter on Edges must be increased.
    The accumulated MFR will be set correct within AccumulateMFR().
*/
void ParMultiGridCL::HandlerTSetPrio( DDD_OBJ obj, DDD_PRIO prio)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    Assert(!(prio==PrioGhost && tp->IsMaster() && !DDD_InfoIsLocal( tp->GetHdr() )),
             DROPSErrCL("HandlerSetPrio: illegal prio for T"),
             DebugParallelC
          );
    if (prio==PrioGhost)
        tp->_Parent= 0;
    if (prio==PrioMaster && tp->GetPrio()==PrioGhost && tp->IsRegularlyRef())
    {   // It may happen, that this routine increases the MFR on a ghost edge, that will be after the transfer
        // not distributed any more. So remeber this tetra and repait the MFR on local edges of this tetra
        // within ParMultiGridCL::AccumulateMFR()
        tp->CommitRegRefMark();
        ToHandleTetra_.push_back(tp);
    }
}


//  Update-Handlers
//-----------------------------
/// \brief Update tetra after recieving
/** link tetra tho faces and set prio to ghost, if no parent is found*/
void ParMultiGridCL::HandlerTUpdate( DDD_OBJ obj)
{
    TetraCL* const tp= ddd_cast<TetraCL*>(obj);
    // link tetra to his faces
    // NOTE: parent has to be linked first, if this Tetra has PrioMaster, has a parent and this parent hasn't linked before
    if (tp->IsMaster() && tp->GetParent() && !tp->GetParent()->GetFace(0)->IsLinkedTo( tp->GetParent() ) )
         // link parent
        for (Uint i=0; i<NumFacesC; ++i)
            const_cast<FaceCL*>(tp->GetParent()->GetFace(i))->LinkTetra( tp->GetParent() );
    if (tp->IsMaster() && !tp->GetParent() && tp->GetLevel()!=0)
    {
        PrioChange(tp, PrioGhost);
        // This case happens within the migration phase. Imagine the following situation:
        //     Proc A has a master copy of a tetra that should be moved to another proc (tetra is
        //         marked for removement on this Proc A)
        //     Proc B moves a ghost copy of the tetra to proc A.
        // Now proc A does not delete the tetrahedron, because it should be recieved by proc B. The both priorities
        // (Ghost and Master) are merged. And Master is wrongly the winner. This is corrected here!
//         AllComment("["<<ProcCL::MyRank()<<"] ====> Set Prio of Tetra " << tp->GetGID() <<" to Ghost, because this tetra has no parent on this proc! (of: This should be OK!)"<<std::endl, ~0);
    }

    // now link this Tetra to his faces
    AllComment("  * Tetra with "<<tp->GetGID()<<" is updated (linked to faces) on proc "<<ProcCL::MyRank()<<std::endl,DebugParallelHardC);
    for (Uint i=0; i<NumFacesC; ++i)
    {
        const_cast<FaceCL*>(tp->GetFace(i))->LinkTetra( tp);
        AllComment("    + to face with GID: " << const_cast<FaceCL*>(tp->GetFace(i))->GetGID()<<std::endl, DebugParallelHardC);
    }
}



/****************************************************************************
* D E L E T E  O B J                                                        *
*****************************************************************************
* If DDD tries to delete Objects there occurs an error, because the STL     *
* does not allow that other code delete objects from a list, vector or so   *
* so we uses function to delete the objects!                                *
****************************************************************************/
/// \brief DDD may use this function to delete a simplex, but this must not be happend!
void ParMultiGridCL::DeleteObj(void * /* buffer*/, size_t /*size*/, int ddd_typ)
{
    std::cerr << "Deleting Object of type " << ddd_typ << " is still missing!" << std::endl;
}



/****************************************************************************
* S E T H A N D L E R A L L                                                 *
*****************************************************************************
* This procedure tells DDD how to treat DROPS-stuff                         *
****************************************************************************/
void ParMultiGridCL::SetAllHandler()
{
    // How to construct Simplices
    DDD_SetHandlerCONSTRUCTOR (VertexCL::GetType(), &HandlerVConstructorC);
    DDD_SetHandlerCONSTRUCTOR (EdgeCL::GetType(),   &HandlerEConstructorC);
    DDD_SetHandlerCONSTRUCTOR (FaceCL::GetType(),   &HandlerFConstructorC);
    DDD_SetHandlerCONSTRUCTOR (TetraCL::GetType(),  &HandlerTConstructorC);

    // How to delete Simplices
    DDD_SetHandlerDELETE(VertexCL::GetType(), &HandlerVDeleteC);
    DDD_SetHandlerDELETE(EdgeCL::GetType(),   &HandlerEDeleteC);
    DDD_SetHandlerDELETE(FaceCL::GetType(),   &HandlerFDeleteC);
    DDD_SetHandlerDELETE(TetraCL::GetType(),  &HandlerTDeleteC);

//  DDD_SetHandlerXFERDELETE(VertexCL::GetType(), &HandlerVDeleteC);
//  DDD_SetHandlerXFERDELETE(EdgeCL::GetType(),   &HandlerEDeleteC);
//  DDD_SetHandlerXFERDELETE(FaceCL::GetType(),   &HandlerFDeleteC);
//  DDD_SetHandlerXFERDELETE(TetraCL::GetType(),  &HandlerTDeleteC);



    // How to transfer simplices
    DDD_SetHandlerXFERCOPY (VertexCL::GetType(), &HandlerVXferC);
    DDD_SetHandlerXFERCOPY (EdgeCL::GetType(),   &HandlerEXferC);
    DDD_SetHandlerXFERCOPY (FaceCL::GetType(),   &HandlerFXferC);
    DDD_SetHandlerXFERCOPY (TetraCL::GetType(),  &HandlerTXferC);

    // How to pack data to a transfer of a simplex
    DDD_SetHandlerXFERGATHER(VertexCL::GetType(), &HandlerVGatherC);
    DDD_SetHandlerXFERGATHER(EdgeCL::GetType(),   &HandlerEGatherC);
    DDD_SetHandlerXFERGATHER(TetraCL::GetType(),  &HandlerTGatherC);

    // How to unpack data from a transfer of a simplex
    DDD_SetHandlerXFERSCATTER(VertexCL::GetType(), &HandlerVScatterC);
    DDD_SetHandlerXFERSCATTER(EdgeCL::GetType(),   &HandlerEScatterC);
    DDD_SetHandlerXFERSCATTER(TetraCL::GetType(),  &HandlerTScatterC);

    // How to handle Tetra right after the transfer!
    DDD_SetHandlerUPDATE     (TetraCL::GetType(), &HandlerTUpdateC);
    DDD_SetHandlerOBJMKCONS  (TetraCL::GetType(), &HandlerTObjMkConsC);
    DDD_SetHandlerSETPRIORITY(TetraCL::GetType(), &HandlerTSetPrioC);

    Comment("- All Handlers are set" << std::endl, DebugParallelC);
}


/****************************************************************************
* G E T  M G                                                                *
****************************************************************************/
/// \brief Recieve a reference to the stored MultiGrid
MultiGridCL& ParMultiGridCL::GetMG()
{
    Assert(_mg!=0, DROPSErrCL("ParMultiGridCL: GetMG: No MultiGrid is assigned"), DebugParallelC);
    return *_mg;
}


/****************************************************************************
* D E B U G I N F O                                                         *
****************************************************************************/
/// \brief Show an simplex by GID
/** Search on proc (or all procs) for a simplex of given GID and show DebugInfo*/
void ParMultiGridCL::Show( DDD_GID gid, char *mesg, int proc)
{
    for (Uint l= 0; l<= _mg->GetLastLevel(); ++l)
    {
        ShowSimplex( _mg->GetVerticesBegin(l), _mg->GetVerticesEnd(l), gid, mesg, proc);
        ShowSimplex( _mg->GetEdgesBegin(l), _mg->GetEdgesEnd(l), gid, mesg, proc);
        ShowSimplex( _mg->GetFacesBegin(l), _mg->GetFacesEnd(l), gid, mesg, proc);
        ShowSimplex( _mg->GetTetrasBegin(l), _mg->GetTetrasEnd(l), gid, mesg, proc);
    }

    if (proc != -1 && proc == ProcCL::MyRank() )
    {
        DDD_HDR hdr= DDD_SearchHdr( gid);
        if (!hdr)
        {
            std::cerr << "...stored only locally." << std::endl;
            return;
        }
        for( int* proclist= DDD_InfoProcList( hdr); *proclist!=-1; proclist+= 2)
            std::cerr << "...stored on proc " << *proclist << " with prio " << PrioToString((Uint)proclist[1]) << std::endl;
    }
}


/// \brief Writes all vertices, edges, faces and tetraeders onto the ostream
void ParMultiGridCL::DebugInfo(std::ostream &os)
{
    os << "I have:\n";
    os << _VertCont->size() << " vertices, " << _EdgeCont->size() << " edges, " << _FaceCont->size() << " faces,"
            << _TetraCont->size() << " tetras" << std::endl;

    os << "\nThe Vertices are:\n";
    MultiGridCL::const_VertexIterator vit=_mg->GetAllVertexBegin();
    for (; vit!=_mg->GetAllVertexEnd(); ++vit){
        vit->DebugInfo(os);
//         if (vit->Unknowns.Exist())
//         {
//             for (size_t i=0; i<_VecDesc.size(); ++i)
//             {
//                 const Uint idx=_VecDesc[i]->RowIdx->GetIdx();
//                 if (vit->Unknowns.Exist(idx))
//                 {
//                     IdxT sysnum=vit->Unknowns(idx);
//                     if (!vit->Unknowns.UnkRecieved(idx))
//                         os << " Unknowns of idx "<<idx<<" at pos "<<sysnum<<" out of Vector: "<< _VecDesc[i]->Data[sysnum]<<std::endl;
//                 }
//                 else
//                     os << " No Unknowns of idx "<<idx<<std::endl;
//             }
//         }
    }

    os << "\nThe Edges are:\n";
    MultiGridCL::const_EdgeIterator eit=_mg->GetAllEdgeBegin();
    for (; eit!=_mg->GetAllEdgeEnd(); ++eit)
    {
        eit->DebugInfo(os);
//         if (eit->Unknowns.Exist())
//         {
//             for (size_t i=0; i<_VecDesc.size(); ++i)
//             {
//                 const Uint idx=_VecDesc[i]->RowIdx->GetIdx();
//                 if (eit->Unknowns.Exist(idx))
//                 {
//                     IdxT sysnum=eit->Unknowns(idx);
//                     if (!eit->Unknowns.UnkRecieved(idx))
//                         os << " Unknowns of idx "<<idx<<" out of Vector: "<< _VecDesc[i]->Data[sysnum]<<std::endl;
//                 }
//                 else
//                     os << " No Unknowns of idx "<<idx<<std::endl;
//             }
//         }
    }

    os << "\nThe Faces are:\n";
    MultiGridCL::const_FaceIterator fit=_mg->GetAllFaceBegin();
    for (; fit!=_mg->GetAllFaceEnd(); ++fit)
        fit->DebugInfo(os);

    os << "\nThe Tetras are:\n";
    MultiGridCL::const_TetraIterator sit(_mg->GetAllTetraBegin());
    for (; sit!=_mg->GetAllTetraEnd(); ++sit)
        sit->DebugInfo(os);

    os << "\n";
}

/// \brief Calc Balance over procs over tetras in the last triangulation level
/** This proc compares the highest number of tetras in the last triangulation level
    with the number of the lowest number of tetras and return the ratio
    (1 is perfect, 0 is unbalanced) */
double ParMultiGridCL::GetBalance()
{
    int myTetras = _mg->_TriangTetra.size();                                    // all procs count their tetras
    int maxTetras = GlobalMax(myTetras);
    int minTetras = GlobalMin(myTetras);

    return std::fabs(1- (double)(maxTetras)/(double)(minTetras));
}

/// \brief Display all types that are defined and declared
void ParMultiGridCL::ShowTypes() const
{
    DDD_TypeDisplay(VertexCL::GetType());
    DDD_TypeDisplay(EdgeCL::GetType());
    DDD_TypeDisplay(FaceCL::GetType());
    DDD_TypeDisplay(TetraCL::GetType());
    DDD_TypeDisplay(AddedScalCL::GetType());
    DDD_TypeDisplay(AddedVecCL::GetType());
    DDD_TypeDisplay(_BndPtT);
    DDD_TypeDisplay(_ChildPtrT);
}

/// \brief Display all interfaces used by ParMultiGridCL
///
/// Show all interfaces, that the DDD-System knows. This procedure uses the DDD-function DDD_IFDisplayAll
void ParMultiGridCL::ShowInterfaces() const
{
    DDD_IFDisplayAll();
}

/// \brief DDD-Consisty-Check
void ParMultiGridCL::ConsCheck()
{
    DDD_ConsCheck();
}

void ParMultiGridCL::ShowTetraIF()
{
    DDD_IFDisplay(_TetraIF);
}

void ParMultiGridCL::ShowEdgeIF()
{
    DDD_IFDisplay(_EdgeIF);
}

/// \brief Display the handled VecDescCL
void ParMultiGridCL::ShowVecDesc()
{
    for (size_t i=0; i<_VecDesc.size(); ++i)
        if (_VecDesc[i])
            std::cout << " - VecDesc["<<i<<"]"
                    << ": Level " << _VecDesc[i]->RowIdx->TriangLevel()
                    << ", #Unknowns " << _VecDesc[i]->RowIdx->NumUnknowns()
                    << ", Idx " <<  _VecDesc[i]->RowIdx->GetIdx()
                    << std::endl;
        else
            std::cout << " - VecDesc["<<i<<"] doen't exists!" << std::endl;
}

/// \brief Get the size of the recieve buffer
size_t ParMultiGridCL::GetRecvBufferSize()
{
    return _RecvBuf.size();
}

/****************************************************************************
* P A R A L L E L   F U N C T I O N S   F R O M   M U L T I G R I D         *
*****************************************************************************
* Declare and Define the simplices and of the Types declared in             *
* parmulitgrid.h                                                            *
****************************************************************************/
void VertexCL::Declare(){
    Assert(!_dddT, DROPSErrCL("VertexCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DDD_TypeDeclare((char*)"Vertex");
}

void EdgeCL::Declare(){
    Assert(!_dddT, DROPSErrCL("EdgeCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DDD_TypeDeclare((char*)"Edge");
}

void FaceCL::Declare(){
    Assert(!_dddT, DROPSErrCL("FaceCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DDD_TypeDeclare((char*)"Face");
}

void TetraCL::Declare(){
    Assert(!_dddT, DROPSErrCL("TetaCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DDD_TypeDeclare((char*)"Tetraeder");
}

void ParMultiGridCL::DeclareBndPtT(){
    Assert(!_BndPtT, DROPSErrCL("ParMultiGridCL: Declare: BndPtT-Type allready declared"), DebugParallelC);
    _BndPtT = DDD_TypeDeclare((char*)"Boundary-Points");
}

void ParMultiGridCL::DeclareChildPtrT(){
    Assert(!_ChildPtrT, DROPSErrCL("ParMultiGridCL: Declare: ChildPtrT-Type allready declared"), DebugParallelC);
    _ChildPtrT = DDD_TypeDeclare((char*)"Tetraeder-Pointer");
}

void AddedScalCL::Declare(){
    Assert(!_dddT, DROPSErrCL("AddedDataCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DDD_TypeDeclare((char*)"Scalar-Unknown-Transfer-Type");
}

void AddedVecCL::Declare(){
    Assert(!_dddT, DROPSErrCL("AddedDataCL: Declare: Type allready declared"), DebugParallelC);
    _dddT = DDD_TypeDeclare((char*)"Vector-Unknown-Transfer-Type");
}


void VertexCL::Define()
{
    VertexCL* v = 0;                                                // example VertexCL for telling DDD some information on this class
    DDD_TypeDefine( _dddT, v,                                       // type and an example
                    EL_LDATA,  &v->_Id, sizeof(v->_Id),
                    EL_GDATA,  &v->_Coord, sizeof(v->_Coord),
                    EL_LDATA,  &v->_BndVerts, sizeof(v->_BndVerts),
                    EL_LDATA,  &v->_RemoveMark, sizeof(v->_RemoveMark),
                    EL_DDDHDR, &v->_dddH,                               // offset header and other members of VertexCL
                    EL_LDATA,  &v->Unknowns, sizeof(v->Unknowns),
                    EL_END,    v+1);
//        EL_DATAPTR, &v->_BndVerts, sizeof(v->_BndVerts),
// WARNUNG: obige Zeile fuehrt dazu, dass der Zeiger _BndVerts ueberschrieben wird; das fuehrt zu Fehlern in HandlerVertexScatter, falls das Zielobjekt schon
//          auf dem Prozess existiert!
}


void EdgeCL::Define()
{
    EdgeCL* e =0;
    DDD_TypeDefine( _dddT, e,
                    EL_DDDHDR, &e->_dddH,
                    EL_OBJPTR, &e->_Vertices, sizeof(e->_Vertices), VertexCL::GetType(),
                    EL_OBJPTR, &e->_MidVertex, sizeof(e->_MidVertex), VertexCL::GetType(),
                    EL_GDATA,  &e->_Bnd, sizeof(e->_Bnd),
                    EL_LDATA,  &e->_MFR, sizeof(e->_MFR),               // Wird mit Interface ... akkumuliert!
                    EL_GDATA,  &e->_AccMFR, sizeof(e->_AccMFR),
                    EL_LDATA,  &e->_RemoveMark, sizeof(e->_RemoveMark),
                    EL_END,    e+1);
    // neglect other members for now:
    // Unknowns
    /* eventuell noch:        EL_GDATA,  &e->_AccMFR, sizeof(e->_AccMFR),*/
}


void FaceCL::Define()
{
    FaceCL* f=0;
    DDD_TypeDefine( _dddT, f,
                    EL_DDDHDR, &f->_dddH,
                    EL_GDATA,  &f->_Bnd, sizeof(f->_Bnd),
                    EL_LDATA,  &f->_RemoveMark, sizeof(f->_RemoveMark),
                    EL_END,    f+1);
    // eventuell noch:        EL_OBJPTR, &f->_Neighbors, sizeof(f->_Neighbors), TetraCL::GetType(),
    // neglect other members for now:
    // Unknowns, _Neighbors
    // _Neighbors sind i.A. auf jedem Proc in anderer Reihenfolge gespeichert!
    //    -> Tetras tragen sich nach Xfer selbst ein.
}


void TetraCL::Define()
{
    TetraCL* t=0;
    DDD_TypeDefine( _dddT, t,
                    EL_DDDHDR, &t->_dddH,
                    EL_GDATA,  &t->_RefRule, sizeof(t->_RefRule),
                    EL_GDATA,  &t->_RefMark, sizeof(t->_RefMark),
                    EL_OBJPTR, &t->_Vertices, sizeof(t->_Vertices), VertexCL::GetType(),
                    EL_OBJPTR, &t->_Edges, sizeof(t->_Edges), EdgeCL::GetType(),
                    EL_OBJPTR, &t->_Faces, sizeof(t->_Faces), FaceCL::GetType(),
                    EL_OBJPTR, &t->_Parent, sizeof(t->_Parent), TetraCL::GetType(),
                    EL_LDATA,  &t->Unknowns, sizeof(t->Unknowns),
                    EL_END,    t+1);
    // neglect other members for now:
    // _Children, Unknowns
}

void ParMultiGridCL::DefineBndPtT()
{
    BndPointCL* b=0;
    DDD_TypeDefine( _BndPtT, b,
                    EL_GDATA,  &b->_BndIdx,  sizeof(b->_BndIdx),        // The Index and
                    EL_GDATA,  &b->_Coord2D, sizeof(b->_Coord2D),       // Coordinate of the Boundary-Point shopuld be submittet
                    EL_END,    b+1);
}

void ParMultiGridCL::DefineChildPtrT()
{
    TetraCL** chp=0;
    DDD_TypeDefine( _ChildPtrT, chp,
                    EL_OBJPTR, &*chp, sizeof(*chp), TetraCL::GetType(),
                    EL_END,    chp+1);
    //delete chp;
}

void AddedScalCL::Define()
{
    AddedScalCL *add=0;
    DDD_TypeDefine( _dddT, add,
                    EL_GDATA, &add->idxVecDesc_, sizeof(add->idxVecDesc_),
                    EL_GDATA, &add->data_, sizeof(add->data_),
                    EL_END,   add+1);
}

void AddedVecCL::Define()
{
    AddedVecCL *add=0;
    DDD_TypeDefine( _dddT, add,
                    EL_GDATA, &add->idxVecDesc_, sizeof(add->idxVecDesc_),
                    EL_GDATA, &add->data_, sizeof(add->data_),
                    EL_END,   add+1);
}

void PrintMG(const DROPS::ParMultiGridCL& pmg, int type)
/** Write out all information about vertices, edges, faces and tetrahedra that
    can be get via the DebugInfo member function of these classes
    \param pmg  the parallel multigrid
    \param type after refinement or migration
*/
{
    const int me=DROPS::ProcCL::MyRank();
    static int REFnum=0;
    static int MIGnum=0;
    char filename[30];

    if (type==REF)
        std::sprintf(filename, "output/%i_MG_REF_%i.mg",me,REFnum++);
    else if (type == MIG)
        std::sprintf(filename, "output/%i_MG_MIG_%i.mg",me,MIGnum++);
    if (me==0)
        std::cout << " - Writing multigrid into: " << filename<< " (for proc 0)"<<std::endl;

    std::ofstream file(filename);
    pmg.DebugInfo(file);
    file.close();
}

bool CheckParMultiGrid(const DROPS::ParMultiGridCL& pmg)
/**  Check parallel and sequential multigrid on each processor
    \param pmg The parallel multigrid
*/
{
    char dat[30];
    std::ostream output (std::cerr.rdbuf());
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream checkfile(dat);
    if (!checkfile){
        IF_MASTER
          std::cerr << "Cannot open file "<<dat<<" to write sanity check output. Using std::cerr"<<std::endl;
    }
    else{
        output.rdbuf(checkfile.rdbuf());
    }
    bool pmg_sane = pmg.IsSane(output),
         mg_sane  = pmg.GetMG().IsSane(output);
    return DROPS::Check(pmg_sane && mg_sane);
}

} // end of namespace DROPS

#endif
