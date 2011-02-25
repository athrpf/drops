/// \file geomselect.tpp
/// \brief offers build/create routines for some standard domains
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Christoph Lehrenfeld, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

namespace DROPS {

template< class T>
void readBoundary ( std::vector<T>& bndcond, const std::string& bnd_type)
{
    size_t idx;
    std::string bndtype( bnd_type), delim("!");

    while ((idx= bndtype.find_first_of( delim)) != std::string::npos )
        bndtype[idx]= ' ';

    std::istringstream bnd( bndtype);

    for( size_t i=0; i<bndcond.size(); ++i){
        bnd >> bndcond[i];
    }
}

template< class BoundaryT>
void BuildBoundaryData( MultiGridCL* &mgp, BoundaryT* &bnddata,
        const std::string& bnd_type_string, const std::string& bnd_funcs_string, match_fun periodic_match)
{
    const BoundaryCL& bnd= mgp->GetBnd();
    const BndIdxT num_bnd= bnd.GetNumBndSeg();

    typedef typename BoundaryT::bnd_val_fun BndValFunT;

    std::vector<int>         bnd_idx  ( num_bnd);
    std::vector<std::string> bnd_names( num_bnd);
    BndCondT* bnd_types = new BndCondT[num_bnd];


    BndValFunT* bnd_fun = new BndValFunT[num_bnd];

    readBoundary( bnd_idx, bnd_type_string);
    readBoundary( bnd_names, bnd_funcs_string);

    typedef std::vector<BoundaryCL::BndType>   BndTypeCont;    
    BndTypeCont perbndt(num_bnd);  

    for( size_t i=0; i<num_bnd; ++i){
        bnd_types[i] = DROPS::BndCondT( bnd_idx[i]);
        bnd_fun[i] = SingletonMapCL<BndValFunT>::getInstance()[bnd_names[i]];
        switch (bnd_types[i])
        {
            case DROPS::Per1BC:
                perbndt[i] = BoundaryCL::Per1Bnd; break;
            case DROPS::Per2BC:
                perbndt[i] = BoundaryCL::Per2Bnd; break;
            default:
                perbndt[i] = BoundaryCL::OtherBnd; break;
        }
    }

    bnddata = new BoundaryT( num_bnd, bnd_types, bnd_fun);
    if (periodic_match)
      bnd.SetPeriodicBnd( perbndt, periodic_match);
    
    delete[] bnd_types;
    delete[] bnd_fun;
}

}// end of namespace DROPS
