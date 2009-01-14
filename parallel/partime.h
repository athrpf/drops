/***************************************************************************
*  File:    partime.h                                                      *
*  Content: This class implements a parallel time-stamping                 *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           10.11.2005                                             *
*  last modified:   10.11.2005                                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file partime.h
/// \brief Parallel time-stamping

#include "parallel/parallel.h"
#include "misc/utils.h"
#include "num/spmat.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

#ifndef _PARTIMERCL_H_
#define _PARTIMERCL_H_
namespace DROPS{

/****************************************************************************
* P A R  T I M E R  C L A S S                                               *
****************************************************************************/
/// \brief Class for parallel time-stamping
/****************************************************************************
* P A R  T I M E R  C L A S S                                               *
****************************************************************************/
class ParTimerCL : public TimerCL
{
  public:
    typedef TimerCL  base_;
  private:
    double   maxduration_;                          // maximal time over the procs
    VectorCL durations_;                            // durations of all procs

    bool calcMax_,                                  // Flag for remembering, if maximal time is calculated
         calcDur_;                                  // and if the array of durations is set

  public:
    ParTimerCL();                                   // constructor

    inline double GetMyTime() const;                            ///< get local time
    inline double GetMaxTime();                                 ///< get time of the proc, that used most time
    inline double GetTime() {return GetMaxTime();}              ///< get time of the proc, that used most time
    inline void   Reset(double time=0);                         ///< reset timer and start measurement
    inline void   Start();                                      ///< start time measurement
    inline void   Stop();                                       ///< stop time measurement
    void PrintAllTime(std::ostream&, int proc=Drops_MasterC);   ///< Proc "proc" print the time of all procs onto the given output stream
};

using std::string;
/****************************************************************************
* T I M E  S T O R E  C L A S S                                             *
****************************************************************************/
/// \brief Class for storing and handling of time measurements
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
/** Perform a barrier and using base class to start time measurement. */
{
    maxduration_=-1.;
    calcMax_=false;
    calcDur_=false;
    ProcCL::Barrier();
    base_::Start();
}

void ParTimerCL::Stop()
/** Use base class to stop time measurement. */
{
    base_::Stop();
    maxduration_=-1.;
    calcMax_=false;
    calcDur_=false;
}

void ParTimerCL::Reset(double time)
/** Reset the time and start the time measurement. */
{
    maxduration_=-1.;
    calcMax_=false;
    calcDur_=false;
    ProcCL::Barrier();
    base_::Reset(time);         // Reset performs also a Start()!
}

double ParTimerCL::GetMyTime() const
{
    return base_::GetTime();
}

double ParTimerCL::GetMaxTime()
/** If maximum over all procs is not computed so far, comput the maximum of times over all procs */
{
    if (!calcMax_){
        maxduration_ = GlobalMax(base_::GetTime());
        calcMax_=true;
    }
    return maxduration_;
}

double TimeStoreCL::AddTime(size_t timeT, double time)
/** Add the given time to the specified time
    \param timeT add to wich time
    \param time  used time
    \return added time
*/
{
    if(timeT>=times_)
        std::cerr << "---> Wrong access on TimeStoreCL!" << std::endl;
    else
        durations_[timeT] += time;
    return time;
}

template<typename TimerT>
double TimeStoreCL::AddTime(const size_t timeT, TimerT timer)
/** Stop the timer and add the time to the specified time
    \param timeT add to wich time
    \param timer  timer
    \return added time
*/
{
    if(timeT>=times_){
        std::cerr << "---> Wrong access on TimeStoreCL!" << std::endl;
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
        std::cerr << "---> Wrong access on TimeStoreCL!" << std::endl;
    else
        return durations_[timeT];
    return 1e300;
}


} // end of namespace DROPS
#endif
