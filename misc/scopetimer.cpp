/// \file scopetimer.cpp
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

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <map>

#include "misc/scopetimer.h"

RealTimerCL::RealTimerCL(double time)
	: _t_begin(timestamp()), _t_end(timestamp()), _time(time)
{
}

void RealTimerCL::Reset(double time)
{
   _time = time;
   _t_begin = timestamp();
}

void RealTimerCL::Start()
{
   _t_begin = timestamp();
}

void RealTimerCL::Stop()
{
   _t_end = timestamp();
   _time += _t_end - _t_begin;
}

double RealTimerCL::GetTime()
{
   return _time;
}

double RealTimerCL::timestamp()
{
#ifdef WIN32
      LARGE_INTEGER time, freq;
      if (QueryPerformanceCounter(&time) == 0)
      {
         DWORD err = GetLastError();
         LPVOID buf;
         FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
                       NULL, err, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                       (LPTSTR) &buf, 0, NULL);
         printf("QueryPerformanceCounter() failed with error %d: %s\n", err, buf);
         exit(1);
      }
      if (QueryPerformanceFrequency(&freq) == 0)
      {
         DWORD err = GetLastError();
         LPVOID buf;
         FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
                       NULL, err, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                       (LPTSTR) &buf, 0, NULL);
         printf("QueryPerformanceFrequency() failed with error %d: %s\n", err, buf);
         exit(1);
      }
      return = Li2Double(time) / Li2Double(freq);
#else
      struct timeval tv;
      gettimeofday(&tv, (struct timezone*)0);
      return ((double)tv.tv_sec + (double)tv.tv_usec / 1000000.0 );
#endif
}



class ScopeTimeCollectorCL
{
private:
   struct SScopeData
   {
      double dTime;
      long   lCall;

      SScopeData& update(double time, long call)
      {
         dTime += time;
         lCall += call;
         return *this;
      }
   };

   typedef std::map<std::string, SScopeData> mapping;

   mapping *_hm_times;
   int _threads;

public:
   ScopeTimeCollectorCL()
   {
#ifdef _OPENMP
      _threads = omp_get_max_threads();
      omp_init_lock(&globalLock);
#else
      _threads = 1;
#endif
      _hm_times = new mapping[_threads];
   }

   ~ScopeTimeCollectorCL()
   {
      mapping::iterator it;
      for (int i=0; i<_threads; i++)
      {
         if (_hm_times[i].size() > 0)
         {
            fprintf(stderr, "\nCollected Timers:\n");
            for(it=_hm_times[i].begin(); it != _hm_times[i].end(); it++)
            {
               fprintf(stderr, "    %38s [%2d](%4ld) :: %12.3f\n",
                   it->first.c_str(), i, it->second.lCall, it->second.dTime);
            }
         }
      }
      delete[] _hm_times;
#ifdef _OPENMP
      omp_destroy_lock(&globalLock);
#endif
   }

   void add( const std::string& name, double time)
   {
#ifdef _OPENMP
      omp_set_lock(&globalLock);
      _hm_times[omp_get_thread_num()][name] = _hm_times[omp_get_thread_num()][name].update(time, (long)1);
      omp_unset_lock(&globalLock);
#else
      _hm_times[0][name] = _hm_times[0][name].update(time, (long)1);
#endif
   }

#ifdef _OPENMP
   omp_lock_t   globalLock;
#endif
};
ScopeTimeCollectorCL scopetimecollector;

ScopeTimerCL::ScopeTimerCL(const std::string& name)
{
   _name = name;
   _timer.Start();
}

ScopeTimerCL::~ScopeTimerCL()
{
   _timer.Stop();
   scopetimecollector.add(_name, _timer.GetTime());
}
