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
