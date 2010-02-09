/// \file pardistributeddatat.h
/// \brief This is the interface between DROPS and a library to handle distributed data, such as DDD
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier, Alin Bastea

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

//How to use: 1. Read the comments before trying to edit this file.
//            2. Comments precede the data they describe.


#include <stdio.h>
#include "geom/simplex.h"
#include "misc/container.h"
#include "parallel/addeddata.h"
#include "parallel/distributeddatatypes.h"


#ifndef DROPS_PARDISTRIBUTEDDATA_H
#define DROPS_PARDISTRIBUTEDDATA_H


namespace DROPS
{

/// \brief Interface for handling dynamic distributed data
/**    Here are declared and defined all the members of the wrapper.
       You can find the implementation of this class (wrapper) in the parddd.cpp file*/
class DynamicDataInterfaceCL
{
    // The declaration of the class' members
  public:
    /// \name Declaring all members
    //@{
    /// \brief Declare a dynamic distributed data type
    static TypeT TypeDeclare (char *name);

    /// \brief Starts transfer phase.
    static void XferBegin();

    /// \brief
    static void XferEnd();

    /// \brief
    static void IFExchange(IFT, size_t, ComProcPtrT, ComProcPtrT);

    /// \brief
    static void IFAOnewayX(IFT, ATTRT, IF_DIRT, size_t, ComProcXPtrT,ComProcXPtrT);

    /// \brief
    static void IFOneway(IFT,IF_DIRT,size_t, ComProcPtrT,ComProcPtrT);

    /// \brief
    static void IFAOneway(IFT ,ATTRT ,IF_DIRT ,size_t , ComProcPtrT,ComProcPtrT);

    /// \brief
    static void IFAExecLocal (IFT, ATTRT, ExecProcPtrT);

    /// \brief
    static void IdentifyObject (HDRT, PROCT, HDRT);

    /// \brief
    static void IFExecLocal (IFT, ExecProcPtrT);

    /// \brief
    static void IdentifyNumber (HDRT, PROCT, int);

    /// \brief
    static IFT IFDefine (int, TypeT *, int, PrioT *, int, PrioT *);

    /// \brief
    static void IFSetName (IFT, char *);

    /// \brief
    static void XferCopyObj (HDRT, PROCT, PrioT);

    /// \brief
    static void HdrDestructor (HDRT);

    /// \brief
    static void XferAddData (int, TypeT);

    /// \brief
    static void SetHandlerCONSTRUCTOR (TypeT, HandlrContrctrT);

    /// \brief
    static void SetHandlerDELETE (TypeT , HandlrDltT);

    /// \brief
    static void SetHandlerXFERCOPY (TypeT, HandlerXfrcpT);

    /// \brief
    static void SetHandlerXFERGATHER (TypeT, HandlerXfrGthrT);

    /// \brief
    static void SetHandlerXFERSCATTER (TypeT, HandlerXfrSctrT);

    /// \brief
    static void SetHandlerUPDATE (TypeT, HandlerUpdtT);

    /// \brief
    static void SetHandlerOBJMKCONS (TypeT, HandlerObjMkConsT);

    /// \brief
    static void SetHandlerSETPRIORITY (TypeT, HandlerSetPrioT);

    /// \brief
    static HDRT SearchHdr (GIDT);

    /// \brief
    static int * InfoProcList (HDRT);

    /// \brief
    static void TypeDisplay (TypeT);

    /// \brief
    static void IFDisplayAll (void);

    /// \brief Returns total #errors
    static int ConsCheck (void);

    /// \brief
    static void IFDisplay (IFT);

    /* This is an elegant implementation of a wrapper around the variable number of
     * arguments function DDD_TypeDefine. However the macros work only on the GCC compiler !!!*/

    /*#define FORCE_INLINE __attribute__((__always_inline__))

    static FORCE_INLINE void TypeDefine (TypeT t,...)
    {
        printf("%d ",__builtin_va_arg_pack_len ());
        DDD_TypeDefine(t, __builtin_va_arg_pack ());
    }*/

    /// \name These functions are the separate initialization of the DDD_TypeDefine method
    ///       They are used to define types for transfering distributed data. Combining all eight 
    ///       functions to a single one is desirable. This is quite ugly ... 
   //@{
    static void TypeDefineAddedVec(DDD_TYPE&, DROPS::AddedVecCL*&, DDD_ELEM_TYPE, DROPS::Uint*, long unsigned int, DDD_ELEM_TYPE, DROPS::Point3DCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::AddedVecCL*);
    static void TypeDefineAddedScal(DDD_TYPE&, DROPS::AddedScalCL*&, DDD_ELEM_TYPE, DROPS::Uint*, long unsigned int, DDD_ELEM_TYPE, double*, long unsigned int, DDD_ELEM_TYPE, DROPS::AddedScalCL*);
    static void TypeDefineChildPtrT(DDD_TYPE&, DROPS::TetraCL**&, DDD_ELEM_TYPE, DROPS::TetraCL**, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::TetraCL**);
    static void TypeDefineBndPtT(DDD_TYPE&, DROPS::BndPointCL*&, DDD_ELEM_TYPE, DROPS::BndIdxT*, long unsigned int, DDD_ELEM_TYPE, DROPS::Point2DCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::BndPointCL*);
    static void TypeDefineTetra(DDD_TYPE&, DROPS::TetraCL*&, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, DROPS::Usint*, long unsigned int, DDD_ELEM_TYPE, DROPS::Usint*, long unsigned int, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::VertexCL*, 4u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::EdgeCL*, 6u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::FaceCL*, 4u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::TetraCL**, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::UnknownHandleCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::TetraCL*);
    static void TypeDefineFace(DDD_TYPE&, DROPS::FaceCL*&, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, const DROPS::BndIdxT*, long unsigned int, DDD_ELEM_TYPE, bool*, long unsigned int, DDD_ELEM_TYPE, DROPS::FaceCL*);
    static void TypeDefineEdge(DDD_TYPE&, DROPS::EdgeCL*&, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::VertexCL*, 2u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::VertexCL**, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::SArrayCL<short unsigned int, 2u>*, long unsigned int, DDD_ELEM_TYPE, short int*, long unsigned int, DDD_ELEM_TYPE, short int*, long unsigned int, DDD_ELEM_TYPE, bool*, long unsigned int, DDD_ELEM_TYPE, DROPS::EdgeCL*);
    static void TypeDefineVertex(DDD_TYPE&, DROPS::VertexCL*&, DDD_ELEM_TYPE, DROPS::IdCL<DROPS::VertexCL>*, long unsigned int, DDD_ELEM_TYPE, DROPS::Point3DCL*, long unsigned int, DDD_ELEM_TYPE, std::vector<DROPS::BndPointCL, std::allocator<DROPS::BndPointCL> >**, long unsigned int, DDD_ELEM_TYPE, bool*, long unsigned int, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, DROPS::UnknownHandleCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::VertexCL*);
    //@}

    /// \brief
    static void SetHandlerXFERDELETE (TypeT, HandlerXferDltT);

    /// \brief
    static PROCT InfoProcs (void);

    /// \brief
    static void Init (int *, char ***);

    /// \brief
    static PROCT InfoMe (void);

    /// \brief Logoff from DDD
    static void Exit (void);

    /// \brief
    static void IdentifyBegin (void);

    /// \brief
    static RETT IdentifyEnd (void);
    //@}

};  // end of the class (the wrapper)

}   // end of namespace
#endif
