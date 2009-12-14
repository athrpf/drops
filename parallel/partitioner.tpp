/***************************************************************************
*  File:    partitioner.h                                                  *
*  Content: - Interface for partitioners                                   *
*           - Metis partitioner                                            *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           18.09.2009                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file partitioner.tpp

namespace DROPS{

template <typename IdxT>
const CSRPartitionerCL<IdxT>::PartitionT& MetisPartitionerCL<IdxT>::PartSerial()
{
    base::part_.resize( base::GetNumVertices());

    int    wgtflag    = 3,                  // Weights on vertices and adjacencies are given
           numflag    = 0,                  // numbering of verts starts by 0 (C-Style)
           nparts     = ProcCL::Size(),     // number of subdomains (per proc one)
           n          = (int)this->GetNumVertices(),
           options[5] = {0,0,0,0,0};        // default options

    if (this->GetMethod()==KWay)
        METIS_PartGraphKway(      &n, xadj_, adjncy_, vwgt_, adjwgt_, &wgtflag, &numflag,  &nparts, options,&edgecut_, part_);
    else if (meth==Recursive)
        METIS_PartGraphRecursive( &n, xadj_, adjncy_, vwgt_, adjwgt_, &wgtflag, &numflag,  &nparts, options,&edgecut_, part_);


    return part_;
}


}   // end of namespace
