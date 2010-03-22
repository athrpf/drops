/// \file fastmarch.tpp
/// \brief fast marching method for reparametrization
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

namespace DROPS{

/** Update the perpendicular foot p to a given dof. If the dof already holds a
    perpendicular foot, the new one is taken, if the the distance is smaller
    \param dof  degree of freedom of the level set unknown
    \param dist new distance
    \param p    perpendicular foot
*/
void ReparamDataCL::UpdatePerp( const IdxT dof, const double dist, const Point3DCL& p)
{
    Assert( UsePerp(), DROPSErrCL("ReparamDataCL::UpdatePerp: Not ready for handling perpendicular feet"), DebugNumericC);
    if ( perpFoot[dof]!=0 && phi.Data[dof]>dist){
        *perpFoot[dof]= p;
    }
    else if ( perpFoot[dof]==0){
        perpFoot[dof]= new Point3DCL(p);
    }
}

/** Compute gradient of level set (as P1 FE on child => const gradient on child)
 */
template <>
inline void InitZeroP1CL<1>::ComputeOnChild( IdxT* Numb, int*, const ChildDataCL& data, LocalP2CL<>& PhiLoc)
{
    SMatrixCL<3,4> p1grad;
    double det;
    Point3DCL pt[4];
    SVectorCL<4> ls;
    for (int v= 0; v < 4; ++v) {
        pt[v]= base::data_.coord[ Numb[data.Vertices[v]]];
        ls[v]= PhiLoc[ data.Vertices[v]];
    }

    P1DiscCL::GetGradients( p1grad, det, pt);
    const Point3DCL n( p1grad*ls);
    const double absdet  = std::abs(det),
                 normgrad= n.norm();

    for (int v= 0; v < 4; ++v) {
        const IdxT MapNr= base::data_.Map( Numb[ data.Vertices[v]]);
#pragma omp critical
{
        sumNormGradPhi_[MapNr]+= normgrad*absdet;
        sumVol_[MapNr]+= absdet;
        base::data_.typ[MapNr]= data_.Finished;
}
    }
}

template<>
inline void InitZeroP1CL<0>::ComputeOnChild(IdxT* Numb, int* sign, const ChildDataCL& data, LocalP2CL<>& PhiLoc)
{
    Point3DCL Schnitt[4];
    int num = 0;
    // Berechnung der Schnittpunkte mit den Kanten des Tetra
    for (int vert = 0; vert < 4; ++vert)
        if (sign[ data.Vertices[vert]] == 0)
            Schnitt[num++] = base::data_.coord[ Numb[data.Vertices[vert]]];

    for (int edge = 0; edge < 6 && num < 4; ++edge) {
        const Ubyte v1 = data.Vertices[ VertOfEdge(edge, 0)],
                    v2 = data.Vertices[ VertOfEdge(edge, 1)];
        if (sign[v1]*sign[v2]==-1) // Vorzeichenwechsel auf edge
        {
            const IdxT Nr1 = Numb[v1], Nr2 = Numb[v2];
            const double bary = InterfacePatchCL::EdgeIntersection(v1, v2, PhiLoc);
            Schnitt[num++] = (1 - bary) * base::data_.coord[Nr1] + bary * base::data_.coord[Nr2];
        }
    }

    if (num < 3)
        throw DROPSErrCL("FastMarchCL::InitZero: intersection missing");

    for (int repeat = 0; repeat < num - 2; ++repeat) { // fuer num==4 (Schnitt ABDC ist viereckig)
        // zwei Dreiecke ABC + DBC betrachten
        if (repeat)
            Schnitt[0] = Schnitt[3];

        const Point3DCL a = Schnitt[1] - Schnitt[0], b = Schnitt[2] - Schnitt[0];

        for ( int vert = 0; vert < 4; ++vert) {
            if ( sign[data.Vertices[vert]] == 0)
                continue;

            const IdxT Nr = Numb[data.Vertices[vert]], MapNr = base::data_.Map(Nr);
            const Point3DCL Crd = base::data_.coord[Nr], c = Crd - Schnitt[0];
            double dist = std::min(c.norm(), (Crd - Schnitt[1]).norm());
            dist = std::min(dist, (Crd - Schnitt[2]).norm());

            double bary1 = inner_prod(a, c) / a.norm_sq(),
                   bary2 = inner_prod(b, c) / b.norm_sq();
            data_.Normalize(bary1, bary2);
            const Point3DCL lotfuss = (1-bary1-bary2) * Schnitt[0] + bary1*Schnitt[1]+bary2*Schnitt[2];
            if (data_.UsePerp()){
                if ( bary1>0 && bary2>0){   // "lotfuss" is located inside tetra
                    if ( data_.perpFoot[MapNr]==0)
                        data_.perpFoot[MapNr]= new Point3DCL( lotfuss);
                    else
                        *data_.perpFoot[MapNr]= lotfuss;
                }
                else{
                    if ( data_.perpFoot[MapNr]!=0) {
                        delete data_.perpFoot[MapNr]; data_.perpFoot[MapNr]=0;
                    }
                }
            }
            dist = std::min(dist, (lotfuss - Crd).norm());

#pragma omp critical
{
            if (base::data_.typ[MapNr] != data_.Finished) {
                base::data_.typ[MapNr] = data_.Finished;
                base::data_.phi.Data[MapNr] = dist;
            }
            else
                base::data_.phi.Data[MapNr] = std::min(dist, base::data_.phi.Data[MapNr]);
}
        }
    }
}

template <int scale>
void InitZeroP1CL<scale>::Perform()
{
#ifdef _PAR
    if (data_.per)
        throw DROPSErrCL("InitZeroP1CL<scale>: Sorry, Periodic boundary conditions are not yet supported by the parallel version");
#endif
    const Uint idx= data_.per ? data_.augmIdx->GetIdx() : data_.phi.RowIdx->GetIdx(),
               lvl= data_.phi.GetLevel();
    int        sign[10];
    int        num_sign[3]; // - 0 +
    IdxT       Numb[10];

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    VecDescCL oldv( data_.phi);
    LocalP2CL<> PhiLoc;

#pragma omp parallel for private( sign, num_sign, Numb, PhiLoc)
    for ( int i=0; i<std::distance(data_.mg.GetTriangTetraBegin(lvl), data_.mg.GetTriangTetraEnd(lvl)); ++i ){
        MultiGridCL::TriangTetraIteratorCL it= data_.mg.GetTriangTetraBegin(lvl)+i;
        for ( int v=0; v<10; ++v){ // collect data on all DoF
            if (v<4)
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            else
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);

            const IdxT MapNr= data_.Map(Numb[v]);
            sign[v]= std::abs( data_.old[MapNr])<1e-8 ? 0 : ( data_.old[MapNr]>0 ? 1 : -1);
            if (sign[v]==0){
#pragma omp critical
                data_.typ[MapNr]= data_.Finished;
            }
        }
        PhiLoc.assign( *it, oldv, NoBndDataCL<>());

        for ( Uint ch=0; ch<MaxChildrenC; ++ch){
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for ( Uint vert= 0; vert<NumVertsC; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];

            const bool intersec= (num_sign[0]*num_sign[2]!=0); // Vorzeichenwechsel

            if (!intersec)
                continue;

            ComputeOnChild( Numb, sign, data, PhiLoc);
        }
    }

    if (scale==1) {
#ifdef _PAR
        // accumulate sums
        if (data_.per) {
            data_.augmIdx->GetEx().Accumulate( sumNormGradPhi_);
            data_.augmIdx->GetEx().Accumulate( sumVol_);
        }
        else {
            data_.phi.RowIdx->GetEx().Accumulate( sumNormGradPhi_);
            data_.phi.RowIdx->GetEx().Accumulate( sumVol_);
        }
#endif
#pragma omp parallel for
        for (int i=0; i<(int)data_.phi.Data.size(); ++i){
            if ( data_.typ[i]== data_.Finished) {                      // dof at interface
                const double scaling= sumNormGradPhi_[i]/sumVol_[i];    // average norm
                data_.phi.Data[i]= std::abs( data_.old[i])/scaling;
            }
        }
    }
}

#ifdef _PAR
/** If level set dof are located on the simplex and distance triangles are associated to
    the corresponding dof, put them into the "tosend"-list for sharing processes
*/
template<typename SimplexT>
int ParInitZeroExactCL::ExecGatherDistTriang(OBJT obj)
{
    SimplexT* const sp= ddd_cast<SimplexT*>(obj);
    const Uint idx= actualData_->phi.RowIdx->GetIdx();
    if ( sp->Unknowns.Exist() && sp->Unknowns.Exist( idx)){
        const IdxT dof= sp->Unknowns( idx);
        if ( !(*distTriangs_)[dof].empty()){
            const ExchangeCL::ProcNumCT neighs= actualData_->phi.RowIdx->GetEx().GetProcs(dof);
            for ( typename ExchangeCL::ProcNumCT::const_iterator it=neighs.begin(); it!=neighs.end(); ++it)
                toSend_->insert( (*distTriangs_)[dof].begin(), (*distTriangs_)[dof].end());
            maxTriangsPerDOF_= std::max( maxTriangsPerDOF_, (*distTriangs_)[dof].size());
        }
    }
    return 0;
}

/** Put send position in the buffer. The first number of the buffer represents the process*/
template<typename SimplexT>
int ParInitZeroExactCL::HandlerDistTriangGatherPos(OBJT objp, void* buf)
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    IdxT* buffer      = static_cast<IdxT*>(buf);
    *(buffer++)= (IdxT) ProcCL::MyRank();
    const Uint idx= actualData_->phi.RowIdx->GetIdx();
    if ( sp->Unknowns.Exist() && sp->Unknowns.Exist( idx)){
        const IdxT dof= sp->Unknowns( idx);
        if ( !(*distTriangs_)[dof].empty()){
            for ( size_t distTriang=0; distTriang<(*distTriangs_)[dof].size(); ++distTriang){
                *(buffer++)= std::distance(toSend_->begin(), toSend_->find((*distTriangs_)[dof][distTriang]));
            }
        }
    }
    *buffer= NoIdx;
    return 0;
}

template<typename SimplexT>
int ParInitZeroExactCL::HandlerDistTriangScatterPos(OBJT objp, void* buf)
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    IdxT*   buffer    = static_cast<IdxT*>(buf);
    const Uint idx= actualData_->phi.RowIdx->GetIdx();
    const int fromProc= *(buffer++);
    if ( sp->Unknowns.Exist() && sp->Unknowns.Exist( idx)){
        const IdxT dof= sp->Unknowns( idx);
        while ( (*buffer)!=NoIdx){
            (*distTriangs_)[dof].push_back( (*actualOffset_)[fromProc]+ *(buffer++));
        }
    }
    return 0;
}
#endif

#ifdef _PAR
template<typename SimplexT>
  int FastmarchingOnMasterCL::HandlerGlobDOFGather(OBJT objp, void* buf)
/** On sender side collect global number of dof on a simplex (if not there send NoIdx)*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    IdxT* buffer      = static_cast<IdxT*>(buf);
    Uint  idx         = actualData_->phi.RowIdx->GetIdx();
    if (sp->Unknowns.Exist() && sp->Unknowns.Exist( idx))
        *buffer= globNumb_[ sp->Unknowns( idx)];      // may be NoIdx if simplex is not exclusive
    return 0;
}

template<typename SimplexT>
  int FastmarchingOnMasterCL::HandlerGlobDOFScatter(OBJT objp, void* buf)
/** On received side collect global number of dof on simplex if sender has send a number*/
{
    SimplexT* const sp = ddd_cast<SimplexT*>(objp);
    const IdxT* buffer = static_cast<IdxT*>(buf);
    Uint  idx          = actualData_->phi.RowIdx->GetIdx();
    if (*buffer!=NoIdx){
        if (sp->Unknowns.Exist() && sp->Unknowns.Exist( idx)){
            Assert(globNumb_[sp->Unknowns( idx)]==NoIdx, DROPSErrCL("FastMarchCL::HandlerGlobDOFScatter: Two exclusive simplices found!"), DebugParallelNumC);
            globNumb_[ sp->Unknowns( idx)]= *buffer;
        }
    }
    return 0;
}

template<typename SimplexT>
  int FastmarchingOnMasterCL::HandlerFinishedGather(OBJT objp, void* buf)
/** On sender-side collect typ and value of distributed DoF*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);              // pointer to simplex
    CoupMarkValST* buffer = static_cast<CoupMarkValST*>(buf);   // pointer to coupling of mark and value
    Uint idx= actualData_->phi.RowIdx->GetIdx();

    if (!sp->Unknowns.Exist() || !sp->Unknowns.Exist( idx))
        return 1;
    Uint Nr= sp->Unknowns( idx);                                // number of DoF
    buffer->mark= actualData_->typ[Nr];                        // set typ
    buffer->val=  actualData_->phi.Data[Nr];                   // set value
    return 0;
}

template<typename SimplexT>
  int FastmarchingOnMasterCL::HandlerFinishedScatter(OBJT objp, void* buf)
/** On receiver side, update value and mark of finished DoFs*/
{
    SimplexT* const sp= ddd_cast<SimplexT*>(objp);
    CoupMarkValST* buffer = static_cast<CoupMarkValST*>(buf);
    Uint idx= actualData_->phi.RowIdx->GetIdx();

    if (!sp->Unknowns.Exist() || !sp->Unknowns.Exist( idx))
        return 1;

    // if a finished dof is transfered
    if ( buffer->mark==actualData_->Finished) {
        Uint Nr= sp->Unknowns( idx);
        // if dof is marked before transfer as finished or within transfer phase
        if ( actualData_->typ[Nr]==actualData_->Finished)
            actualData_->phi.Data[Nr] = std::min( actualData_->phi.Data[Nr], buffer->val);
        else
            actualData_->phi.Data[Nr] = buffer->val;
        // set dof as finished
        actualData_->typ[Nr]= actualData_->Finished;
    }

    return 0;
}
#endif

#ifdef _PAR
template<class SimplexT>
  int ParDirectDistanceCL::HandlerFrontierGather(OBJT objp, void* buffer)
/** Put value of phi, the perpendicular foot (if available) and proc id
    into the buffer. If the simplex is not marked as frontier, take -1
    as proc id
*/
{
    SimplexT* sp   = ddd_cast<SimplexT*>( objp);                // simplex (vertex/edge)
    TransferST *buf= static_cast<TransferST*>( buffer);         // content of message
    IdxT dof= sp->Unknowns( actualData_->phi.RowIdx->GetIdx());   // where to find phi
    // fill buffer
    buf->value = actualData_->phi.Data[dof];
    buf->perp  = actualData_->perpFoot[dof] ? *(actualData_->perpFoot[dof]) : Point3DCL(std::numeric_limits<double>::max());
    buf->procID= actualData_->typ[dof]==ReparamDataCL::Finished ? ProcCL::MyRank() : -1;
    return 0;
}

template<class SimplexT>
  int ParDirectDistanceCL::HandlerFrontierScatter(OBJT objp, void* buffer)
{
    SimplexT* sp      = ddd_cast<SimplexT*>( objp);
    TransferST* buf   = static_cast<TransferST*>( buffer);
    IdxT dof= sp->Unknowns( actualData_->phi.RowIdx->GetIdx());
    if ( buf->procID>=0)
        onProc_[dof].push_front( TransferST(*buf));
    return 0;
}

#endif

InitZeroP2CL::RepTetra::RepTetra( const TetraCL& t, LocalP1CL<Point3DCL>* Gref, const ReparamDataCL& data) 
    : w2b(t)
{
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( G, Gref, T);
    const Uint idx= data.phi.RowIdx->GetIdx();
    for ( Uint v=0; v<NumVertsC; ++v){
        const IdxT dof= t.GetVertex(v)->Unknowns(idx);
        coord[v] = data.coord[dof];
        valPhi[v]= data.phi.Data[dof];
    }
    for ( Uint e=0; e<NumEdgesC; ++e){
        const IdxT dof= t.GetEdge(e)->Unknowns(idx);
        coord[e+NumVertsC] = data.coord[dof];
        valPhi[e+NumVertsC]= data.phi.Data[dof];
    }
    baryCenter= GetBaryCenter( t);
}

}   // end of namespace
