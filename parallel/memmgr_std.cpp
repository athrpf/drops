/// \file memmgr_std.cpp
/// \brief memory management for DDD
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
/// this is a modificated file from DDD

/****************************************************************************/
/*                                                                          */
/* File:      memmgr_std.c                                                  */
/*                                                                          */
/* Purpose:   basic memory management module                                */
/*            (with standard malloc() calls)                                */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/04/27 kb  begin                                            */
/*            96/01/20 kb  updated to DDD V1.5                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#include <cstdlib>
#include "ppif.h"

#ifdef __GNUC__
#  define __UNUSED__ __attribute__((__unused__))
#else
#  define __UNUSED__
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


extern "C"
{

void *memmgr_AllocPMEM (std::size_t size)
{
    return std::malloc(size);
}


void memmgr_FreePMEM (void *buffer)
{
    std::free(buffer);
}


void *memmgr_AllocOMEM (std::size_t size, __UNUSED__ int ddd_typ,
                        __UNUSED__ int proc, __UNUSED__ int attr)
{
    return std::malloc(size);
}

/*extern void deleteObjC(void *buffer, size_t size, int ddd_typ);*/

void memmgr_FreeOMEM (__UNUSED__ void *buffer, __UNUSED__ std::size_t size,
                      __UNUSED__ int ddd_typ)
{
    /*deleteObjC(buffer, size, ddd_typ);*/
    /*#ifndef _PAR*/
    /*free(buffer);*/
    /*#endif*/
}


void *memmgr_AllocAMEM (std::size_t size)
{
    return std::malloc(size);
}


void memmgr_FreeAMEM (void *buffer)
{
    std::free(buffer);
}


void *memmgr_AllocTMEM (std::size_t size)
{
    return std::malloc(size);
}


void memmgr_FreeTMEM (void *buffer)
{
    std::free(buffer);
}


void *memmgr_AllocHMEM (std::size_t size)
{
    return std::malloc(size);
}


void memmgr_FreeHMEM (void *buffer)
{
    std::free(buffer);
}


void memmgr_MarkHMEM (void)
{
}


void memmgr_ReleaseHMEM (void)
{
}


void memmgr_Init (void)
{
}

}
