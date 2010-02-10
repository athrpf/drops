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
       You can find the implementation of this class (wrapper) in the parddd.cpp file
       */
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

    /// \brief End of transfer phase.
    static void XferEnd();

    /// \brief Exchanges data via an interface
    static void IFExchange(IFT, size_t, ComProcPtrT, ComProcPtrT);

    static void IFAOnewayX(IFT, ATTRT, IF_DIRT, size_t, ComProcXPtrT,ComProcXPtrT);

    /// \brief Transfer one way
    static void IFOneway(IFT,IF_DIRT,size_t, ComProcPtrT,ComProcPtrT);

    static void IFAOneway(IFT ,ATTRT ,IF_DIRT ,size_t , ComProcPtrT,ComProcPtrT);

    static void IFAExecLocal (IFT, ATTRT, ExecProcPtrT);

    /// \brief Indentifying object which are created by different processes
    static void IdentifyObject (HDRT, PROCT, HDRT);

    static void IFExecLocal (IFT, ExecProcPtrT);

    static void IdentifyNumber (HDRT, PROCT, int);

    static IFT IFDefine (int, TypeT *, int, PrioT *, int, PrioT *);

    /// \brief Allows to define a textual description for interfaces
    static void IFSetName (IFT, char *);

    /// \brief Transfer-command for copying a local DDD object to another processor.
    static void XferCopyObj (HDRT, PROCT, PrioT);

    /// \brief Remove object's header from DDD management
    static void HdrDestructor (HDRT);

    /// \brief Transfer array of additional data objects with a DDD local object
    static void XferAddData (int, TypeT);

    /// \brief Creates C++-Objects
    static void SetHandlerCONSTRUCTOR (TypeT, HandlrContrctrT);

    static void SetHandlerDELETE (TypeT , HandlrDltT);

    static void SetHandlerXFERCOPY (TypeT, HandlerXfrcpT);

    static void SetHandlerXFERGATHER (TypeT, HandlerXfrGthrT);

    static void SetHandlerXFERSCATTER (TypeT, HandlerXfrSctrT);

    static void SetHandlerUPDATE (TypeT, HandlerUpdtT);

    static void SetHandlerOBJMKCONS (TypeT, HandlerObjMkConsT);

    static void SetHandlerSETPRIORITY (TypeT, HandlerSetPrioT);

    /// \brief Searching object which are created by different processes
    static HDRT SearchHdr (GIDT);

    static int * InfoProcList (HDRT);

    /// \brief Show defined TypeT
    static void TypeDisplay (TypeT);

    /// \brief Display overview of all DDD interfaces.
    static void IFDisplayAll (void);

    /// \brief Check DDD runtime consistency. Returns total number of errors
    static int ConsCheck (void);

    /// \brief Display overview of single DDD interface.
    static void IFDisplay (IFT);

    /* This is an elegant implementation of a wrapper around the variable number of
     * arguments function DDD_TypeDefine. However the macros work only on the GCC compiler !!!*/

    /*#define FORCE_INLINE __attribute__((__always_inline__))

    static FORCE_INLINE void TypeDefine (TypeT t,...)
    {
        printf("%d ",__builtin_va_arg_pack_len ());
        DDD_TypeDefine(t, __builtin_va_arg_pack ());
    }*/

    /// \brief these are the callings of the variable number of arguements function DDD_TypeDefine
    static void TypeDefineAddedVec(DDD_TYPE&, DROPS::AddedVecCL*&, DDD_ELEM_TYPE, DROPS::Uint*, long unsigned int, DDD_ELEM_TYPE, DROPS::Point3DCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::AddedVecCL*);

    static void TypeDefineAddedScal(DDD_TYPE&, DROPS::AddedScalCL*&, DDD_ELEM_TYPE, DROPS::Uint*, long unsigned int, DDD_ELEM_TYPE, double*, long unsigned int, DDD_ELEM_TYPE, DROPS::AddedScalCL*);

    static void TypeDefineChildPtrT(DDD_TYPE&, DROPS::TetraCL**&, DDD_ELEM_TYPE, DROPS::TetraCL**, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::TetraCL**);

    static void TypeDefineBndPtT(DDD_TYPE&, DROPS::BndPointCL*&, DDD_ELEM_TYPE, DROPS::BndIdxT*, long unsigned int, DDD_ELEM_TYPE, DROPS::Point2DCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::BndPointCL*);

    static void TypeDefineTetra(DDD_TYPE&, DROPS::TetraCL*&, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, DROPS::Usint*, long unsigned int, DDD_ELEM_TYPE, DROPS::Usint*, long unsigned int, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::VertexCL*, 4u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::EdgeCL*, 6u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::FaceCL*, 4u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::TetraCL**, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::UnknownHandleCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::TetraCL*);

    static void TypeDefineFace(DDD_TYPE&, DROPS::FaceCL*&, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, const DROPS::BndIdxT*, long unsigned int, DDD_ELEM_TYPE, bool*, long unsigned int, DDD_ELEM_TYPE, DROPS::FaceCL*);

    static void TypeDefineEdge(DDD_TYPE&, DROPS::EdgeCL*&, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, DROPS::SArrayCL<DROPS::VertexCL*, 2u>*, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::VertexCL**, long unsigned int, DDD_TYPE, DDD_ELEM_TYPE, DROPS::SArrayCL<short unsigned int, 2u>*, long unsigned int, DDD_ELEM_TYPE, short int*, long unsigned int, DDD_ELEM_TYPE, short int*, long unsigned int, DDD_ELEM_TYPE, bool*, long unsigned int, DDD_ELEM_TYPE, DROPS::EdgeCL*);

    static void TypeDefineVertex(DDD_TYPE&, DROPS::VertexCL*&, DDD_ELEM_TYPE, DROPS::IdCL<DROPS::VertexCL>*, long unsigned int, DDD_ELEM_TYPE, DROPS::Point3DCL*, long unsigned int, DDD_ELEM_TYPE, std::vector<DROPS::BndPointCL, std::allocator<DROPS::BndPointCL> >**, long unsigned int, DDD_ELEM_TYPE, bool*, long unsigned int, DDD_ELEM_TYPE, DDD_HEADER*, DDD_ELEM_TYPE, DROPS::UnknownHandleCL*, long unsigned int, DDD_ELEM_TYPE, DROPS::VertexCL*);

    static void SetHandlerXFERDELETE (TypeT, HandlerXferDltT);

    static PROCT InfoProcs (void);

    /// \brief Initialisation of the library.
    static void Init (int *, char ***);

    static PROCT InfoMe (void);

    /// \brief Clean-up of the DDD library.
    static void Exit (void);

    /// \brief Begin identification phase.
    static void IdentifyBegin (void);

    /// \brief End of identification phase.
    static RETT IdentifyEnd (void);
    //@}

};  // end of the class (the wrapper)

}   // end of namespace
#endif
