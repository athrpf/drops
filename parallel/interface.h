//**************************************************************************
// File:    interface.h                                                    *
// Content: Classes that handle DDD-Interfaces                             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:    September, 19th 2007                                           *
//          - AllSimplexIFCL added                                         *
// Begin:   Februar, 8th 2006                                              *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file interface.h
/// \brief Classes for handling DDD-Interfaces
#ifndef DROPS_INTERFACE_H
#define DROPS_INTERFACE_H

#include <ddd.h>
#include <iostream>
#include "geom/multigrid.h"

namespace DROPS{

/****************************************************************************
* I N T E R F A C E  C L A S S                                              *
****************************************************************************/
/// \brief Manage the interfaces over simplices that are able to store unknowns
/** In order to access distributed simplices and send messages between
    distributed simplices, DDD provides interfaces. This class manages
    interfaces bewteen simplices that may store unknowns.
*/
/****************************************************************************
* I N T E R F A C E  C L A S S                                              *
****************************************************************************/
template <class SimplexT>
class InterfaceCL
{
  private:
    static DDD_IF SimplexIF_;

  public:
    InterfaceCL() {}
    ~InterfaceCL() {}

    static void InitIF();                   // Init Interface over Master-Simplices
    static inline DDD_IF GetIF();           // Get Interface over Simplices
    static void ShowIF();                   // Display Simplex-Interface
};



/****************************************************************************
* A L L  I N T E R F A C E  C L A S S                                       *
****************************************************************************/
/// \brief Manage the interfaces over simplices
/** In order to access distributed simplices and send messages between
    distributed simplices, DDD provides interfaces. This class manages
    interfaces bewteen all distributed simplices
*/
/****************************************************************************
* A L L  I N T E R F A C E  C L A S S                                       *
****************************************************************************/
template<typename SimplexT>
class AllSimplexIFCL
{
  private:
    static DDD_IF AllSimplexIF_;

  public:
    AllSimplexIFCL() {}
    ~AllSimplexIFCL() {}

    static void InitIF();                   // Init interface
    static DDD_IF GetIF();                  // Get Interface
    static void ShowIF();                   // Display interface
};


// init of static members
//------------------------
template <class SimplexT> DDD_IF InterfaceCL<SimplexT>::SimplexIF_=0;
template <class SimplexT> DDD_IF AllSimplexIFCL<SimplexT>::AllSimplexIF_=0;


// template and inline functions
//--------------------------------

/// \brief Get the Interface over Simplices
template <class SimplexT>
  DDD_IF InterfaceCL<SimplexT>::GetIF()
{
    return SimplexIF_;
}


/// \brief Define the interfaces just over simplices that are able to store unknowns
template <class SimplexT>
  void InterfaceCL<SimplexT>::InitIF()
{
    if (SimplexIF_!=0){
        std::cout << "=====> InterfaceCL: Allready defined that Interface!" << std::endl;
    }

    DDD_TYPE  O[1];
    DDD_PRIO  A[1], B[1];

    O[0] = SimplexT::GetType();
    A[0] = PrioHasUnk;// A[1]= PrioVGhost; A[2]= PrioGhost;
    B[0] = PrioHasUnk;// B[1]= PrioVGhost; B[2]= PrioGhost;

    SimplexIF_ = DDD_IFDefine(1, O, 1, A, 1, B); // only Master
    char name[40];
    std::sprintf(name, "Assemble IF for Type %d", SimplexT::GetType());
    DDD_IFSetName( SimplexIF_, name);
}


/// \brief Print information about interface onto std output
template <class SimplexT>
  void InterfaceCL<SimplexT>::ShowIF()
{
    DDD_IFDisplay(SimplexIF_);
}


/// \brief Get the Interface over Simplices
template <class SimplexT>
  DDD_IF AllSimplexIFCL<SimplexT>::GetIF()
{
    Assert(AllSimplexIF_!=0, DROPSErrCL("AllSimplexIFCL::GetIF: Interface not init"), DebugParallelNumC);
    return AllSimplexIF_;
}


/// \brief Init Interface over all distributed types of simplices
template <class SimplexT>
  void AllSimplexIFCL<SimplexT>::InitIF()
{
    if (AllSimplexIF_!=0){
        std::cout << "=====> InterfaceCL: Allready defined that Interface!" << std::endl;
    }

    DDD_TYPE  O[1];
    DDD_PRIO  A[4], B[4];

    O[0] = SimplexT::GetType();
    A[0] = PrioHasUnk; A[1]= PrioMaster; A[2]= PrioGhost; A[3]=PrioVGhost;
    B[0] = PrioHasUnk; B[1]= PrioMaster; B[2]= PrioGhost; B[3]=PrioVGhost;

    AllSimplexIF_ = DDD_IFDefine(1, O, 4, A, 4, B);
    char name[50];
    std::sprintf(name, "All Simplex Interface for type %d", SimplexT::GetType());
    DDD_IFSetName( AllSimplexIF_, name);
}


/// \brief Print information about interface onto std output
template <class SimplexT>
  void AllSimplexIFCL<SimplexT>::ShowIF()
{
    DDD_IFDisplay(AllSimplexIF_);
}


}   // end of namespace DROPS
#endif
