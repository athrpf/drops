//**************************************************************************
// File:    utils.h                                                        *
// Content: useful stuff that fits nowhere else                            *
// Author:  Joerg Peters, Volker Reichelt, IGPM RWTH Aachen                *
// Version: 0.8                                                            *
// History: begin - September, 16, 2000                                    *
//                                                                         *
// Remarks:                                                                *
//**************************************************************************

#ifndef DROPS_UTILS_H
#define DROPS_UTILS_H


// #include <limits> // TODO: Do we have limits with gcc-snapshots or SGI-CC?
#include <functional>
#include <algorithm>
#include <iosfwd>
#include <string>
#include <ctime>


namespace DROPS
{


typedef unsigned int      Uint;
typedef unsigned long int Ulint;
typedef unsigned short    Usint;
typedef signed char       byte;
typedef unsigned char     Ubyte;


const double DoubleEpsC = 1.0e-9; // numeric_limits<double>::epsilon();

// constants that switch on debugging for specified portions of code
#define DebugContainerC      1
#define DebugRefineEasyC     2
#define DebugRefineHardC     4
#define DebugNumericC        8
#define DebugUnknownsC      16
#define DebugNoReuseSparseC 32

#define cdebug std::cout

//#define DROPSDebugC 25  //(DROPS::DebugNumericC | DROPS::DebugUnknownsC | DROPS::DebugContainerC )
//#define DROPSDebugC ~0  // all bits set
#define DROPSDebugC 0 

#if DROPSDebugC
#  define Assert(a,b,c) (_Assert((a),(b),(c)))
#else
#  define Assert(a,b,c) ((void)0)
#endif

#if DROPSDebugC
#  define Comment(a,b) do { if ((b) & DROPSDebugC) cdebug << a; } while (false)
#else
#  define Comment(a,b) ((void)0)
#endif

#ifdef __GNUC__
#  define __UNUSED__ __attribute__((__unused__))
#else
#  define __UNUSED__
#endif


#undef DROPS_EXP_TEMPL_CONSTR_FOR_VALARRAY_DERIVATIVE
#define DROPS_EXP_TEMPL_CONSTR_FOR_VALARRAY_DERIVATIVE(theClass, thebase_type) \
template <class X__>                                                           \
  explicit theClass (const X__& x__): thebase_type( x__) {}

#undef  DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE
#define DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE(theOp, theClass, theT, thebase_type) \
theClass& operator theOp (const theT s)                                                  \
{                                                                                        \
    *static_cast<thebase_type*>( this) theOp s; return *this;                            \
}                                                                                        \
template <class VT>                                                                      \
  inline theClass<theT>&                                                                 \
  operator theOp (const VT& v)                                                           \
{                                                                                        \
    Assert( this->size()==v.size(), #theClass #theOp ": incompatible dimensions",        \
        DebugNumericC);                                                                  \
    *static_cast<thebase_type*>( this) theOp v; return *this;                            \
}


#undef  DROPS_ASSIGNMENT_OPS_FOR_VALARRAY_DERIVATIVE
#define DROPS_ASSIGNMENT_OPS_FOR_VALARRAY_DERIVATIVE(theClass, theT, thebase_type) \
/*assignment*/                                                                     \
DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE(=, theClass, theT, thebase_type)       \
/*computed assignment*/                                                            \
DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE(+=, theClass, theT, thebase_type)      \
DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE(-=, theClass, theT, thebase_type)      \
DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE(*=, theClass, theT, thebase_type)      \
DROPS_ASSIGNMENT_OP_FOR_VALARRAY_DERIVATIVE(/=, theClass, theT, thebase_type)      \

#undef  DROPS_DEFINE_VALARRAY_DERIVATIVE
#define DROPS_DEFINE_VALARRAY_DERIVATIVE(theClass, theT, thebase_type)     \
/*The expression template constructor*/                                    \
DROPS_EXP_TEMPL_CONSTR_FOR_VALARRAY_DERIVATIVE(theClass, thebase_type)     \
/*assignment and computed assignment*/                                     \
DROPS_ASSIGNMENT_OPS_FOR_VALARRAY_DERIVATIVE(theClass, theT, thebase_type)


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
_Assert(A assertion, E exc, Uint DebugLevel=~0)
{
    if (DebugLevel&DROPSDebugC) 
        if (!assertion) throw exc;
}


template <class A>
inline void
_Assert(A assertion, const char* msg, Uint DebugLevel=~0)
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

    static void ResetCounter(Ulint i= 0) { _Counter= i; }

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


// This is needed, as the C++-standard comittee deemed selectors for
// pairs unneccessary.
// Function-object that returns the second component of a pair.
template <class Pair>
struct select2nd : public std::unary_function<Pair, typename Pair::second_type>
{
  typename Pair::second_type& operator()(Pair& x) const {
    return x.second;
  }
  const typename Pair::second_type& operator()(const Pair& x) const {
    return x.second;
  }
};


template <class T>
class ptr_iter
{
  public:
    typedef std::random_access_iterator_tag iterator_category;
    typedef T                               value_type;
    typedef ptrdiff_t                       difference_type;
    typedef T*                              pointer;
    typedef T&                              reference;

  private:
    T** p_;

  public:
    // default copy-ctor, copy-assignment-op and dtor
    ptr_iter (T** p) : p_( p) {}
    
    reference operator*  () const { return **p_; }
    pointer   operator-> () const { return *p_; }

    ptr_iter        operator+ (difference_type d) const { return p_ + d; }
    difference_type operator- (ptr_iter b)        const { return p_ - b.p_; }

    reference operator[] (difference_type d) const { return **(p_ + d); }

    ptr_iter& operator++ ()    { ++p_; return *this; }
    ptr_iter& operator++ (int) { ptr_iter tmp( p_); ++p_; return tmp; }
    ptr_iter& operator-- ()    { --p_; return *this; }
    ptr_iter& operator-- (int) { ptr_iter tmp( p_); --p_; return tmp; }

    friend bool operator== (const ptr_iter& a, const ptr_iter& b) { return a.p_ == b.p_; }
    friend bool operator!= (const ptr_iter& a, const ptr_iter& b) { return a.p_ != b.p_; }
    friend bool operator<  (const ptr_iter& a, const ptr_iter& b) { return a.p_ < b.p_; }
};


template <class T>
struct ConstHelperCL
{
    typedef       T stripped_type;
    typedef const T const_type;
};

template <class T>
struct ConstHelperCL<const T>
{
    typedef       T stripped_type;
    typedef const T const_type;
};


} // end of namespace DROPS

#endif
