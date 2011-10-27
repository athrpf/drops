/// \file scopetimer.h
/// \brief scope timer
/// \author LNM RWTH Aachen: SC RWTH Aachen: Christian Terboven
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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/
#ifndef SCOPETIMER_H
#define SCOPETIMER_H

#ifdef WIN32
 #include <Windows.h>
 #undef min
 #undef max
 #define Li2Double(x) ((double)((x).HighPart) * 4.294967296E9 + (double)((x).LowPart))
#else
 #include <sys/time.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include <map>
#include <string>



class RealTimerCL
{
private:
   double _t_begin, _t_end;
   double _time;

public:
   RealTimerCL(double time=0.);
   void Reset(double time= 0);
   void Start();
   void Stop();
   double GetTime();
   
   double timestamp();
};



//*******************************************************************
// S c o p e T i m e r C L
//   accumulates elapsed times between creation and 
//   destruction of objects with the same 'name'
//*******************************************************************
class ScopeTimerCL
{
private: 
   std::string _name;
   RealTimerCL _timer;

public:
   ScopeTimerCL(const std::string& name);
   ~ScopeTimerCL();
};
typedef ScopeTimerCL ScopeTimer;



#endif
