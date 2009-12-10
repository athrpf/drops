/// \file partime.cpp
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
#include "parallel/partime.h"

namespace DROPS{

ParTimerCL::ParTimerCL() : durations_(ProcCL::Size())
{
    maxduration_=-1.;                                                                       // Flag for the maximal time
    calcMax_=false;                                                                         // maxtime
    calcDur_=false;                                                                         // and the time of all proces are not computed yet
    Start();
}

void ParTimerCL::PrintAllTime(std::ostream &os, int proc)
{
    double myTime = GetMyTime();                                                            // calculate the time this proc needed
    durations_[ProcCL::MyRank()] = myTime;                                                  // and store this time in the array of durations

    if (!calcDur_)
        ProcCL::Gather(myTime, Addr(durations_),  proc);                                    // Proc "proc" collects all times

    if (ProcCL::MyRank()==proc)                                                             // only "proc" should print the results onto "os"
        os << "ParTimerCL gets following results:" << std::endl;

    GetMaxTime();                                                                           // Compute global maximum used time

    for (int i=0; i<ProcCL::Size(); ++i)
        if (ProcCL::MyRank()==proc)
            os << "Proc " << i << ": " << durations_[i] << " sec" << '\n';                  // and print the time of proc "i" onto "os"

    if (ProcCL::MyRank()==proc)
    {
        os << "MaxTime is: " << maxduration_ << " sec" << std::endl;                        // print the maximal time also onto "os"
        calcDur_=true;                                                                      // only proc "proc" has the time of all procs
        calcMax_=true;                                                                      // and the maximal time
    }
}

double ParTimerCL::TestBandwidth(std::ostream& os, int messageSize, int numTest)
/** Measure the bandwidth between processor 0 and processor size-1.
    \param messageSize Size of a single message in MB
    \param numTest Number of messages to be send between two processors
    \return bandwidth in GB/s
 */
{
    if (ProcCL::IamMaster())
        os << "Testing bandwidth ...\n";

    if (ProcCL::Size()==1){
        if (ProcCL::IamMaster())
            os << " => Skipping, because only a single processor is involved" << std::endl;
        return 0.;
    }

    std::valarray<byte> buffer('a', messageSize*1024*1024);
    ParTimerCL timer;
    timer.Start();
    ProcCL::RequestT req;
    if (ProcCL::IamMaster()){
        for (int i=0; i<numTest; ++i){
            req= ProcCL::Isend(buffer, ProcCL::Size()-1, 27);
            ProcCL::Wait(req);
        }
    }
    else if (ProcCL::MyRank()==ProcCL::Size()-1){
        for (int i=0; i<numTest; ++i){
            req= ProcCL::Irecv(buffer, ProcCL::Master(), 27);
            ProcCL::Wait(req);
        }
    }
    timer.Stop();
    const double bandwidth= (messageSize/1024.)*numTest/timer.GetTime();
    if (ProcCL::IamMaster())
        os << " => Measured bandwidth " << bandwidth << " GB/s" << std::endl;

    return bandwidth;
}

TimeStoreCL::TimeStoreCL(size_t times) :
    durations_(times,0), describer_(times," "), times_(times),
    overall_(0.), counter_(0), counter_descr_(" ")
/** Constructs a TimeStorCL, that is able to store a specified number of times
    \param times number of stored times
*/
{
}


void TimeStoreCL::SetDescriber(size_t timeT, string str)
{
    if(timeT>=times_)
        std::cout << "---> Wrong access on TimeStoreCL!" << std::endl;
    else
        describer_[timeT] = str;
}

void TimeStoreCL::Print(std::ostream &os)
{
    size_t max_length=counter_descr_.length();
    max_length = (7>max_length ? 7 : max_length);
    for (size_t i=0; i<times_; ++i)
        max_length = (describer_[i].length()>max_length ? describer_[i].length() : max_length);
    max_length+=5;
    if (ProcCL::IamMaster())
    {
        os.precision(6);
        os << "Times are used for:\n";
        for (size_t i=0; i<times_; ++i)
        {
            os.fill('.'); os.setf(std::ios::fixed);
            os << std::setiosflags(std::ios::left)
               << std::setw(max_length) << describer_[i]
               << '\t' << durations_[i] << '\n';
        }
        os.fill('.'); os.setf(std::ios::fixed);
        os << std::setiosflags(std::ios::left) << std::setw(max_length)
           << "Overall" << '\t' << overall_ << '\n';
        os << std::setiosflags(std::ios::left) << std::setw(max_length)
           << counter_descr_ << '\t' << counter_ << std::endl;
        os.unsetf(std::ios::fixed);
    }
}

} // end of namespace DROPS
