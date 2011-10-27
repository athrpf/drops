/// \file interfacePatch.tpp
/// \brief Computes 2D patches and 3D cuts of tetrahedra and interface
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky; SC RWTH Aachen:

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

namespace DROPS
{
inline double
InterfacePatchCL::EdgeIntersection (Uint v0, Uint v1, const LocalP2CL<> & philoc)
{
    if (LinearEdgeIntersection) return philoc[v0]/(philoc[v0]-philoc[v1]);
    else {
        const double l0= philoc[v0], l1= philoc[v1],
            lm= philoc( AllEdgeBaryCenter_[v0][v1]); // Value of Phi in the edge-barycenter.
        // If l0*l1<0 the quadratic equation with p(0)=l0, p(1)=l1, p(1/2)=lm has exactly one root in (0,1).
        const double quadcoeff= 2.*l0 + 2.*l1 - 4.*lm, lincoeff= 4.*lm - 3.*l0 - l1;
        if ( std::fabs( quadcoeff) < std::fabs( lincoeff)*8.*std::numeric_limits<double>::epsilon()) // linear LS-function
            return l0/(l0 - l1);
        const double rt= std::sqrt( std::pow(4.*lm - (l0 + l1), 2) - 4.*l0*l1);
        const double x0= (-lincoeff - rt)/(2.*quadcoeff),
                     x1= (-lincoeff + rt)/(2.*quadcoeff);
        Assert( (0 < x0 && x0 < 1.) || (0 < x1 && x1 < 1.),
            "InterfacePatchCL::EdgeIntersection: Excessive roundoff-error with quadratic level-set-function",
            DebugNumericC);
        return (0. < x0 && x0 < 1.) ? x0 : x1;
    }
}

template<class ValueT>
ValueT InterfaceTetraCL::quad( const LocalP2CL<ValueT>& f, double absdet, bool part /*bool debug*/)
{
    const ChildDataCL& data= GetChildData( ch_);
    typedef BaryCoordCL* BaryPtrT;
    BaryPtrT BaryPtr[4];
    if (intersec_<3)
    { // cuts = Kind + leere Menge
//if (debug) std::cout <<"Fall <3:\tKind + leere Menge\n";
        if ( part == (num_sign_[2]>0) )
        { // integriere ueber Kind
            for (int i=0; i<4; ++i)
                BaryPtr[i]= &BaryDoF_[data.Vertices[i]];
            return P1DiscCL::Quad( f, BaryPtr)*(ch_ < 8 ? 0.125 : 1.)*absdet;
        }
        else
            return ValueT();
    }
    else if (intersec_==3)
    { // cuts = Tetra + Tetra-Stumpf
//if (debug) std::cout <<"Fall 3:\tTetra + Tetra-Stumpf\n";
        int vertA= -1;  // cut-Tetra = APQR
        const int signA= num_sign_[0]==1 ? -1 : 1; // sign of vert A occurs only once
        for (int i=0; vertA==-1 && i<4; ++i)
            if (sign_[data.Vertices[i]]==signA) vertA= i;
        for (int i=0; i<3; ++i)
            BaryPtr[i]= &Bary_[i];
        BaryPtr[3]= &BaryDoF_[data.Vertices[vertA]];

        const double volFrac= VolFrac( BaryPtr);
        const ValueT quadTetra= P1DiscCL::Quad( f, BaryPtr)*(absdet*volFrac);
//if (debug) std::cout << "vertA = " << vertA << "\tvolFrac = " << volFrac << "\t= 1 / " << 1/volFrac << std::endl;
        if ( part == (signA==1) )
            return quadTetra;
        else // Gesamt-Tetra - cut-Tetra
        {
            for (int i=0; i<4; ++i)
                BaryPtr[i]= &BaryDoF_[data.Vertices[i]];
            return P1DiscCL::Quad( f, BaryPtr)*(ch_ < 8 ? 0.125 : 1.)*absdet - quadTetra;
        }
    }
    else // intersec_==4
    { // cuts = 2x 5-Flaechner mit 6 Knoten. connectivity:     R---------S
      //                                                      /|        /|
      //                                                     A-|-------B |
      //                                                      \|        \|
      //                                                       P---------Q
//if (debug) std::cout <<"Fall 4: 2x 5-Flaechner\t";
        int vertAB[2], // cut mit VZ==part = ABPQRS
            signAB= part ? 1 : -1;
            for (int i=0, k=0; i<4 && k<2; ++i)
            if (sign_[data.Vertices[i]]==signAB) vertAB[k++]= i;
        // connectivity AP automatisch erfuellt, check for connectivity AR
        const bool AR= vertAB[0]==VertOfEdge(Edge_[2],0) || vertAB[0]==VertOfEdge(Edge_[2],1);
//if (debug) if (!AR) std::cout << "\nAR not connected!\n";

//if (debug) std::cout << "vertA = " << vertAB[0] << "\tvertB = " << vertAB[1] << std::endl;
//if (debug) std::cout << "PQRS on edges\t"; for (int i=0; i<4; ++i) std::cout << Edge_[i] << "\t"; std::cout << std::endl;
        // Integriere ueber Tetras ABPR, QBPR, QBSR    (bzw. mit vertauschten Rollen von Q/R)
        // ABPR    (bzw. ABPQ)
        BaryPtr[0]= &BaryDoF_[data.Vertices[vertAB[0]]];
        BaryPtr[1]= &BaryDoF_[data.Vertices[vertAB[1]]];
        BaryPtr[2]= &Bary_[0];    BaryPtr[3]= &Bary_[AR ? 2 : 1];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cout << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        ValueT integral= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        // QBPR    (bzw. RBPQ)
        BaryPtr[0]= &Bary_[AR ? 1 : 2];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cout << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        integral+= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        // QBSR    (bzw. RBSQ)
        BaryPtr[2]= &Bary_[3];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cout << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        integral+= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        return absdet*integral;
    }
}

template<class ValueT>
void InterfaceTetraCL::quadBothParts( ValueT& int_pos, ValueT& int_neg, const LocalP2CL<ValueT>& f, double absdet)
{
    const ChildDataCL& data= GetChildData( ch_);
    typedef BaryCoordCL* BaryPtrT;
    BaryPtrT BaryPtr[4], chTetra[4];

    for (int i=0; i<4; ++i)
        chTetra[i]= &BaryDoF_[data.Vertices[i]];
    const ValueT quadChild= P1DiscCL::Quad( f, chTetra)*(ch_ < 8 ? 0.125 : 1.)*absdet;

    if (intersec_<3)
    { // cuts = Kind + leere Menge
//if (debug) std::cout <<"Fall <3:\tKind + leere Menge\n";
        if ( num_sign_[2]>0 )
        {
            int_pos= quadChild;    int_neg= ValueT();
        }
        else
        {
            int_neg= quadChild;    int_pos= ValueT();
        }
    }
    else if (intersec_==3)
    { // cuts = Tetra + Tetra-Stumpf
//if (debug) std::cout <<"Fall 3:\tTetra + Tetra-Stumpf\n";
        int vertA= -1;  // cut-Tetra = APQR
        const int signA= num_sign_[0]==1 ? -1 : 1; // sign of vert A occurs only once
        for (int i=0; vertA==-1 && i<4; ++i)
            if (sign_[data.Vertices[i]]==signA) vertA= i;
        for (int i=0; i<3; ++i)
            BaryPtr[i]= &Bary_[i];
        BaryPtr[3]= &BaryDoF_[data.Vertices[vertA]];

        const double volFrac= VolFrac( BaryPtr);
        const ValueT quadTetra= P1DiscCL::Quad( f, BaryPtr)*(absdet*volFrac);
//if (debug) std::cout << "vertA = " << vertA << "\tvolFrac = " << volFrac << "\t= 1 / " << 1/volFrac << std::endl;
        if ( signA==1 )
        {
            int_pos= quadTetra;   int_neg= quadChild - quadTetra;
        }
        else
        {
            int_neg= quadTetra;   int_pos= quadChild - quadTetra;
        }
    }
    else // intersec_==4
    { // cuts = 2x 5-Flaechner mit 6 Knoten. connectivity:     R---------S
      //                                                      /|        /|
      //                                                     A-|-------B |
      //                                                      \|        \|
      //                                                       P---------Q
//if (debug) std::cout <<"Fall 4: 2x 5-Flaechner\t";
        int vertAB[2], // cut mit VZ == + = ABPQRS
            signAB= 1;
        for (int i=0, k=0; i<4 && k<2; ++i)
            if (sign_[data.Vertices[i]]==signAB) vertAB[k++]= i;
        // connectivity AP automatisch erfuellt, check for connectivity AR
        const bool AR= vertAB[0]==VertOfEdge(Edge_[2],0) || vertAB[0]==VertOfEdge(Edge_[2],1);
//if (debug) if (!AR) std::cout << "\nAR not connected!\n";
 //std::cout << "vertA = " << vertAB[0] << "\tvertB = " << vertAB[1] << std::endl;
 //std::cout << "PQRS on edges\t"; for (int i=0; i<4; ++i) std::cout << Edge_[i] << "\t"; std::cout << std::endl;
        // Integriere ueber Tetras ABPR, QBPR, QBSR    (bzw. mit vertauschten Rollen von Q/R)
        // ABPR    (bzw. ABPQ)
        BaryPtr[0]= &BaryDoF_[data.Vertices[vertAB[0]]];
        BaryPtr[1]= &BaryDoF_[data.Vertices[vertAB[1]]];
        BaryPtr[2]= &Bary_[0];    BaryPtr[3]= &Bary_[AR ? 2 : 1];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cout << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        ValueT integral= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        // QBPR    (bzw. RBPQ)
        BaryPtr[0]= &Bary_[AR ? 1 : 2];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cout << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        integral+= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        // QBSR    (bzw. RBSQ)
        BaryPtr[2]= &Bary_[3];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cout << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        integral+= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        int_pos= absdet*integral;    int_neg= quadChild - int_pos;
    }
}

template<class ValueT>
ValueT InterfaceTriangleCL::quad2D( const LocalP2CL<ValueT>& f, Uint tri) const
{
    const BaryCoordCL* const Vert= &GetBary(tri); // array of the 3 vertices of triangle i in barycentric coords
    BaryCoordCL Bary;                             // barycenter in barycentric coords
    double sum= 0;

    for (int i=0; i<3; ++i)
    {
        Bary+= Vert[i];
        sum+= f(Vert[i]);
    }
    Bary/= 3.;

    return (0.375*f(Bary) + sum/24.)*GetAbsDet( tri);
}

template <class It>
  std::pair<double, double>
  h_interface (It begin, It end, const VecDescCL& ls)
{
    double hmin= 1e99, hmax= -1.0; //, hmean= 0.;
    // size_t num= 0;
    // EdgeCL* emax= 0;

    for (; begin != end; ++begin)
        if (   (ls.Data[begin->GetVertex( 0)->Unknowns( ls.RowIdx->GetIdx())]
                *ls.Data[begin->Unknowns( ls.RowIdx->GetIdx())] <= 0.)
            || (ls.Data[begin->GetVertex( 1)->Unknowns( ls.RowIdx->GetIdx())]
                *ls.Data[begin->Unknowns( ls.RowIdx->GetIdx())] <= 0.)) {
            const double h= (begin->GetVertex( 0)->GetCoord() - begin->GetVertex( 1)->GetCoord()).norm();
            hmin= std::min( hmin, h);
            // if (h > hmax) emax= &*begin;
            hmax= std::max( hmax, h);
            // hmean+= h;
            // ++num;
        }
    // std::cout << "mean : " << hmean/num << '\n';
    // emax->DebugInfo( std::cout);
    // emax->GetVertex( 0)->DebugInfo( std::cout);
    // emax->GetVertex( 1)->DebugInfo( std::cout);

#ifdef _PAR
    hmin = ProcCL::GlobalMin( hmin);
    hmax = ProcCL::GlobalMax( hmax);
#endif
    return std::make_pair( hmin, hmax);
}

} // end of namespace DROPS
