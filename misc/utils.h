//**************************************************************************
// File:    utils.h                                                        *
// Content: useful stuff that fits nowhere else                            *
// Author:  Joerg Peters, Volker Reichelt, IGPM RWTH Aachen                *
// Version: 0.8                                                            *
// History: begin - September, 16, 2000                                    *
//                                                                         *
// Remarks:                                                                *
//**************************************************************************

#ifndef _UTILS_H_
#define _UTILS_H_


// #include <limits> // TODO: Do we have limits with gcc-snapshots or SGI-CC?
#include <functional>
#include <algorithm>
#include <iosfwd>
#include <string>
#include <time.h>


namespace DROPS
{


typedef unsigned int      Uint;
typedef unsigned long int Ulint;
typedef unsigned short    Usint;
typedef signed char       byte;
typedef unsigned char     Ubyte;


const double DoubleEpsC = 1.0e-9; // numeric_limits<double>::epsilon();

// constants that switch on debugging for specified portions of code
const Uint DebugContainerC  = 1;
const Uint DebugRefineEasyC = 2;
const Uint DebugRefineHardC = 4;
const Uint DebugNumericC    = 8;
const Uint DebugUnknownsC   =16;

#define cdebug std::cout

//#define DROPSDebugC 25  //(DROPS::DebugNumericC | DROPS::DebugUnknownsC | DROPS::DebugContainerC )
//#define DROPSDebugC -1u
#define DROPSDebugC 0 

#if DROPSDebugC
#  define Assert(a,b,c) (_Assert((a),(b),(c)))
#else
#  define Assert(a,b,c) ((void)0)
#endif

#if DROPSDebugC
#  define Comment(a,b) {if ((b) & DROPSDebugC) cdebug << a;}
#else
#  define Comment(a,b) ((void)0)
#endif

template <class In, class T>
inline bool is_in( In beg, In end, const T& value )
{
    return std::find(beg,end,value) != end;
}


template <class In, class Pred>
inline bool is_in_if( In beg, In end, Pred p )
{
    return std::find_if(beg,end,p) != end;
}


template <class In, class Out, class Pred>
Out copy_if( In beg, In end, Out dest, Pred p )
{
    while (beg!=end) { if (p(*beg)) *dest++= *beg; ++beg; }
    return dest;
}


template <class T>
class ref_to_ptr : public std::unary_function<T&, T*>
{
  public:
    T* operator() (T& arg) const { return static_cast<T*>(&arg); }
};


//**************************************************************************
// Class:    DROPSErrCL                                                    *
// Purpose:  base class for all classes that are thrown as exception.      *
// Remarks:  Classes should derive their own (hopefully more powerful)     *
//           error-class and donate meaningful error-messages.             *
//           The default error handler only prints the error message...    *
//**************************************************************************

class DROPSErrCL
{
  protected:
    std::string _ErrMesg;

  public:
    DROPSErrCL() : _ErrMesg("") {}
    DROPSErrCL(const std::string& mesg) : _ErrMesg(mesg) {}
    virtual ~DROPSErrCL() {}

    virtual std::ostream& what  (std::ostream&) const;
    virtual void          handle()              const;
};


template <class E, class A>
inline void
_Assert(A assertion, E exc, Uint DebugLevel=-1)
{
    if (DebugLevel&DROPSDebugC) 
        if (!assertion) throw exc;
}


template <class A>
inline void
_Assert(A assertion, const char* msg, Uint DebugLevel=-1)
{
    if (DebugLevel&DROPSDebugC) 
        if (!assertion) throw DROPSErrCL(msg);
}


//**************************************************************************
// Class:    IdCL                                                          *
// Purpose:  provides a unique identifier for an object.                   *
// Remarks:  We use the template argument to specify the class whose       *
//           objects will carry an Id.                                     *
//**************************************************************************

template <class type>
class IdCL
{
private:
    static Ulint _Counter;

    Ulint _Identity;

public:
    IdCL () : _Identity(_Counter++) {}
    IdCL (Ulint Identity) : _Identity(Identity) {}
    // Default Copy-ctor

    Ulint GetCounter () const { return _Counter; }
    Ulint GetIdent   () const { return _Identity; }

    bool operator == (const IdCL<type>& Id) const
        { return Id._Identity == _Identity; }
    bool operator != (const IdCL<type>& Id) const { return !(*this==Id); }
    bool operator <  (const IdCL<type>& Id) const
        { return _Identity < Id._Identity; }
};

// Initialize the counter only once!
template <class type> Ulint IdCL<type>::_Counter = 0;


//**************************************************************************
//                                                                         *
//     T i m e r C L :   get to know how fast DROPS is !  :-)              *
//                                                                         *
//**************************************************************************

class TimerCL
{
  private:
    clock_t _t_begin, _t_end;
    double _time;          // local timer
    static clock_t _gt_begin, _gt_end;
    static double _gtime;  // global timer
  public:
    TimerCL(double time=0.) : _t_begin(clock()), _t_end(clock()), _time(time) {}
    void Reset(double time= 0) { _time= time; _t_begin= clock(); }
    void Start()               { _t_begin= clock(); }
    void Stop()                { _t_end= clock(); _time+= double(_t_end - _t_begin)/CLOCKS_PER_SEC; }
    double GetTime() const     { return _time; }

    void GReset(double time= 0) { _gtime= time; _gt_begin= clock(); }
    void GStart()               { _gt_begin= clock(); }
    void GStop()                { _gt_end= clock(); _gtime+= double(_gt_end - _gt_begin)/CLOCKS_PER_SEC; }
    double GetGTime() const     { return _gtime; }
};

} // end of namespace DROPS

#endif
