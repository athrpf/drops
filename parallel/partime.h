/// \file partime.h
/// \brief parallel time-stamping
/// \author LNM RWTH Aachen: SC RWTH Aachen: Oliver Fortmeier

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

#include "parallel/parallel.h"
#include "misc/utils.h"
#include "num/spmat.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

#ifndef DROPS_PARTIMER_H
#define DROPS_PARTIMER_H

namespace DROPS{

/****************************************************************************
* P A R  T I M E R  C L A S S                                               *
****************************************************************************/
/// \brief Class for parallel time-stamping
/****************************************************************************
* P A R  T I M E R  C L A S S                                               *
****************************************************************************/
class ParTimerCL
{
  private:
    double   maxduration_;                          // maximal time over the procs
    VectorCL durations_;                            // durations of all procs
    double   start_;
    double   end_;

    bool calcMax_,                                  // Flag for remembering, if maximal time is calculated
         calcDur_;                                  // and if the array of durations is set

  public:
      ParTimerCL();                                   // constructor

      inline double GetMyTime() const;                            ///< get local time
      inline double GetMaxTime();                                 ///< get time of the proc, that used most time
      inline double GetTime() {return GetMaxTime();}              ///< get time of the proc, that used most time
      inline void   Reset();                                      ///< reset timer and start measurement
      inline void   Start();                                      ///< start time measurement
      inline void   Stop();                                       ///< stop time measurement
      void PrintAllTime(std::ostream&, int proc=Drops_MasterC);   ///< Proc "proc" print the time of all procs onto the given output stream
      /// \brief Measure bandwidth
      static double TestBandwidth(std::ostream&, int messageSize=200, int numTest=5);
};

using std::string;
/****************************************************************************
* T I M E  S T O R E  C L A S S                                             *
****************************************************************************/
/// \brief Class for storing and handling of time measurements
/** This class uses MPI timer to compute times */
/****************************************************************************
* T I M E  S T O R E  C L A S S                                             *
****************************************************************************/
class TimeStoreCL
{
  private:
    std::vector<double> durations_;                 // stored times
    std::vector<string> describer_;                 // describer for times
    size_t times_;                                  // number of stored times
    double overall_;                                // overall stored time
    int counter_;                                   // counting something
    string counter_descr_;                          // describer of the counter

  public:
    /// \brief Std-constructor
    TimeStoreCL(size_t times);

    /// \brief Add a time to stored time
    inline double AddTime(size_t timeT, double time);

    /// \brief Add a time to stored time
    template<typename TimerT>
    inline double AddTime(const size_t timeT, TimerT timer);

    /// \brief Set overall time
    void SetOverall(double time) {overall_=time;}

    /// \brief Get overall time
    double GetOverall() const { return overall_; }

    /// \brief Get a stored time
    inline double GetTime(size_t timeT) const;

    /// \brief Increment the counter
    void IncCounter(int i=1) {counter_+=i;}
    /// \brief Get counter
    int  GetCounter() const {return counter_;}

    /// \brief Describe stored time
    void SetDescriber(size_t timeT, string);
    /// \brief Describe counter
    void SetCounterDescriber(string str) {counter_descr_=str;}

    /// \brief Reset all times and describers
    void reset() {
        durations_.resize(0); durations_.resize(times_);
        describer_.resize(0); describer_.resize(times_);
    }

    /// \brief Print information onto an output stream
    void Print(std::ostream&);
};

/****************************************************************************
*   I M P L E M E N T A T I O N   O F   I N L I N E   F U N C T I O N S     *
****************************************************************************/
void ParTimerCL::Start()
/** Perform a barrier and use wall time to start time measurement. */
{
    maxduration_=-1.;
    calcMax_=false;
    calcDur_=false;
    ProcCL::Barrier();
    start_= ProcCL::Wtime();
}

void ParTimerCL::Stop()
/** Stop time measurement. */
{
    end_= ProcCL::Wtime();
    maxduration_=-1.;
    calcMax_=false;
    calcDur_=false;
}

void ParTimerCL::Reset()
/** Reset the time and start the time measurement. */
{
    this->Start();
}

double ParTimerCL::GetMyTime() const
{
      return end_-start_;
}

double ParTimerCL::GetMaxTime()
/** If maximum over all procs has not been computed so far, compute the maximum of times over all procs */
{
    if (!calcMax_){
        maxduration_ = ProcCL::GlobalMax(end_-start_);
        calcMax_=true;
    }
    return maxduration_;
}

double TimeStoreCL::AddTime(size_t timeT, double time)
/** Add the given time to the specified time
    \param timeT add to which time
    \param time  used time
    \return added time
*/
{
    if(timeT>=times_)
        std::cout << "---> Wrong access on TimeStoreCL!" << std::endl;
    else
        durations_[timeT] += time;
    return time;
}

template<typename TimerT>
double TimeStoreCL::AddTime(const size_t timeT, TimerT timer)
/** Stop the timer and add the time to the specified time
    \param timeT add to which time
    \param timer  timer
    \return added time
*/
{
    if(timeT>=times_){
        std::cout << "---> Wrong access on TimeStoreCL!" << std::endl;
        return -1.;
    }
    timer.Stop();
    const double duration= timer.GetTime();
    durations_[timeT] += duration;
    return duration;
}


double TimeStoreCL::GetTime(size_t timeT) const
/** Ask for a stored time
    \param timeT type of stored time
*/
{
    if(timeT>=times_)
        std::cout << "---> Wrong access on TimeStoreCL!" << std::endl;
    else
        return durations_[timeT];
    return 1e300;
}


} // end of namespace DROPS
#endif
