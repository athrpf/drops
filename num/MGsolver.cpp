/// \file
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Maxim Larin, Patrick Esser, IGPM

#include "num/MGsolver.h"

namespace DROPS
{

void MGDataCL::RemoveCoarseResetFinest()
{
    if (this->empty()) return;
    //RemoveCoarse
    MGDataCL::iterator it=this->end();
    --it;
    this->erase(this->begin(), it);
    //ResetFinest
    this->begin()->A.Data.clear();
    this->begin()->P.Data.clear();
    this->begin()->B.Data.clear();
    this->begin()->BT.Data.clear();
    this->begin()->Mpr.Data.clear();
    this->begin()->Mvel.Data.clear();
    this->begin()->PPr.Data.clear();
    this->begin()->AN.Data.clear();
    this->begin()->ABlock = &this->begin()->A.Data;
}

void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end)
{
    Uint lvl= 0;
    for(const_MGDataIterCL coarse= begin, fine=++begin; fine!=end; ++coarse, ++fine, ++lvl)
    {
//        if (lvl==1) std::cerr << fine->A.Data << std::endl;
        const Uint nu= coarse->Idx.NumUnknowns;
        VectorCL ei(nu), Ai(nu);

        std::cerr << nu << " unknowns on level " << lvl << std::endl;


        if (nu>700)
        {
            std::cerr << "Check skipped: too many unknowns, too much time..." << std::endl;
            continue;
        }


        for(Uint i=0; i<nu; ++i)
        {
            if (i!=0) ei[i-1]=0; else ei[nu-1]=0; ei[i]=1;
            Ai= transp_mul( fine->P.Data, fine->A.Data * (fine->P.Data * ei) );
            for(Uint j=0; j<nu; ++j)
            {
                if (std::fabs(Ai[j] - coarse->A.Data(i,j))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": A_H(" << i << ", " << j <<")= "
                              << coarse->A.Data(i,j) << " != " << Ai[j] << std::endl;
                }
            }
        }
    }
}

void CheckStokesMGData( const_MGDataIterCL begin, const_MGDataIterCL end )
{
    Uint lvl= 0;
    for(const_MGDataIterCL coarse= begin, fine=++begin; fine!=end; ++coarse, ++fine, ++lvl)
    {
        const Uint n= coarse->Idx.NumUnknowns;
        const Uint m= coarse->IdxPr.NumUnknowns;
        const Uint nun= n+m;
        VectorCL eui(n), epi(m), Aui(n), Api(m);

        std::cerr << nun << " unknowns on level " << lvl << std::endl;

        if (nun>11000)
        {
            std::cerr << "Check skipped: too many unknowns, too much time..." << std::endl;
            continue;
        }

        std::cerr << "Check velocity" << std::endl;
        eui=0;epi=0;
        for(Uint i=0; i<n; ++i)
        {
            if (i!=0) eui[i-1]=0; else eui[n-1]=0; eui[i]=1;
            Aui= transp_mul( fine->P.Data, VectorCL(fine->A.Data * (fine->P.Data * eui) + transp_mul( fine->B.Data, VectorCL(fine->PPr.Data * epi) )));
            Api= transp_mul( fine->PPr.Data, VectorCL(fine->B.Data * (fine->P.Data * eui )));
//            std::cerr << "    i = " << i << std::endl;
//            std::cerr << "    velocity part " << std::endl;
            for(Uint j=0; j<n; ++j)
            {
                if (std::fabs(Aui[j] - coarse->A.Data(i,j))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": Au_H(" << i << ", " << j <<")= "
                              << coarse->A.Data(i,j) << " != " << Aui[j] << std::endl;
                }
            }
//            std::cerr << "    pressure part" << std::endl;
            for(Uint j=0; j<m; ++j)
            {
                if (std::fabs(Api[j] - coarse->B.Data(j,i))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": Au_H(" << i << ", " << j <<")= "
                              << coarse->B.Data(j,i) << " != " << Api[j] << std::endl;
                }
            }
        }

        std::cerr << "Check pressure" << std::endl;
	eui=0;epi=0;
        for(Uint i=0; i<m; ++i)
        {
            if (i!=0) epi[i-1]=0; else epi[m-1]=0; epi[i]=1;
            Aui= transp_mul( fine->P.Data, VectorCL(fine->A.Data * (fine->P.Data * eui) + transp_mul( fine->B.Data, VectorCL(fine->PPr.Data * epi ))));
            Api= transp_mul( fine->PPr.Data, VectorCL(fine->B.Data * (fine->P.Data * eui )));
//            std::cerr << "    i = " << i << std::endl;
//            std::cerr << "    velocity part " << std::endl;
            for(Uint j=0; j<n; ++j)
            {
                if (std::fabs(Aui[j] - coarse->B.Data(i,j))>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": Au_H(" << i << ", " << j <<")= "
                              << coarse->B.Data(i,j) << " != " << Aui[j] << std::endl;
                }
            }
//            std::cerr << "    pressure part" << std::endl;
            for(Uint j=0; j<m; ++j)
            {
                if (std::fabs(Api[j])>1e-6)
                {
                    std::cerr << "Auf coarse-Level " << lvl << ": Au_H(" << i << ", " << j <<")= "
                              << 0.0 << " != " << Api[j] << std::endl;
                }
            }
        }

        std::cerr << "Check P1 * one_coarse = one_fine" << std::endl; 

        DROPS::VectorCL cones( 1.0, m );
        std::cout << cones <<  std::endl;

        const Uint mf= fine->IdxPr.NumUnknowns;
        DROPS::VectorCL fones( 1.0, mf );
        std::cout << fones <<  std::endl;

        DROPS::VectorCL vtmp = fine->PPr.Data * cones;
        std::cout << vtmp <<  std::endl;

        DROPS::VectorCL vdelta (fones - fine->PPr.Data * cones);
        std::cout << vdelta <<  std::endl;

        std::cout << " vdelta.norm() = " << norm(vdelta) << std::endl;
    }
}


} // end of namespace DROPS
