/// \file distributeddatatypes.h
/// \brief  Redefinintion of the types from the DDD library
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

#include <stdio.h>

#ifndef DROPS_HEADER_H
#define DROPS_HEADER_H

#define USE_DDD

#ifdef _PAR
# ifdef USE_DDD
#   include <ddd.h>
# endif
# ifdef OWN_DDD
#  include <myownddd.h>
# endif
#endif


namespace DROPS
{
#ifdef USE_DDD
    typedef DDD_HEADER HEADERT;                  ///< Object header
    typedef DDD_RET RETT;                        ///< Return types for DDD functions
    typedef DDD_IF IFT;                          ///< Interface type
    typedef ComProcPtr ComProcPtrT;              ///< Function type
    typedef ComProcXPtr ComProcXPtrT;            ///< Function type
    typedef DDD_TYPE TypeT;                      ///< Type of a distributed object
    typedef DDD_ATTR ATTRT;                      ///< attribute of a distributed object
    typedef DDD_IF_DIR IF_DIRT;                  ///< direction of an interface
    typedef ExecProcPtr ExecProcPtrT;            ///< Function type
    typedef DDD_HDR HDRT;                        ///< each object stores a header which can be found here
    typedef DDD_GID GIDT;                        ///< global id type
    typedef DDD_PROC PROCT;                      ///< process numbers
    typedef DDD_PRIO PrioT;                      ///< priority
    typedef DDD_OBJ OBJT;                        ///< pointer to an object
    typedef HandlerCONSTRUCTOR HandlrContrctrT;  ///< handler how to construct an object
    typedef HandlerDELETE HandlrDltT;            ///< Function type
    typedef HandlerXFERCOPY HandlerXfrcpT;       ///< Function type
    typedef HandlerXFERGATHER HandlerXfrGthrT;   ///< Function type
    typedef HandlerXFERSCATTER HandlerXfrSctrT;  ///< Function type
    typedef HandlerUPDATE HandlerUpdtT;          ///< Function type
    typedef HandlerOBJMKCONS HandlerObjMkConsT;  ///< Function type
    typedef HandlerSETPRIORITY HandlerSetPrioT;  ///< Function type
    typedef HandlerXFERDELETE HandlerXferDltT;   ///< Function type
#define MarkHdrInvalid(DDD_HDR) DDD_MarkHdrInvalid(DDD_HDR); ///< Macro definition of the creator of the header
#endif


/// \brief This class is the wrapper around the methods declared in geom/simplex.h
///        which because of the TypeDeclare method had to be included inside parallel/pardistributeddata.h
class DynamicDataInterfaceExtraCL
{
  public:
    /// \brief
    static int InfoIsLocal (HDRT);
    /// \brief
    static int * InfoProcList (DDD_HDR);
    /// \brief
    static void HdrConstructor (DDD_HDR, DDD_TYPE, DDD_PRIO, DDD_ATTR);
    /// \brief
    static void HdrConstructorMove (DDD_HDR, DDD_HDR);
    /// \brief
    static void XferDeleteObj (DDD_HDR);
};
}
#endif
