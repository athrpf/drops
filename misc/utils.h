/// \file
/// \brief Useful stuff that fits nowhere else.

#ifndef DROPS_UTILS_H
#define DROPS_UTILS_H


// #include <limits> ///< \todo Do we have limits with gcc-snapshots or SGI-CC?
#include <functional>
#include <algorithm>
#include <iosfwd>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>


namespace DROPS
{
/// \name Basic types
/// These abbreviations for compound-typenames are used everywhere in DROPS.
//@{
typedef unsigned int      Uint;
typedef unsigned long int Ulint;
typedef unsigned short    Usint;

/// \brief Mainly for tables in topo.h and topo.cpp that store topological data
/// for the refinement algorithm.
typedef signed char       byte;
/// \brief Mainly for tables in topo.h and topo.cpp that store topological data
/// for the refinement algorithm.
typedef unsigned char     Ubyte;
//@}


/// Used in equality-tests for floating point numbers.
const double DoubleEpsC = 1.0e-9; // numeric_limits<double>::epsilon();

/// \name Code-groups for debugging.
/// \brief Constants that group the code for debugging and error-reporting.
//@{
#define DebugContainerC      1
#define DebugRefineEasyC     2
#define DebugRefineHardC     4
#define DebugNumericC        8
#define DebugUnknownsC      16
#define DebugNoReuseSparseC 32
//@}

/// The stream for dedug output.
#define cdebug std::cout

/// \brief This macro controls, for which portions of the code debugging and
/// error-reporting is activated.
//#define DROPSDebugC 25  //(DROPS::DebugNumericC | DROPS::DebugUnknownsC | DROPS::DebugContainerC )
//#define DROPSDebugC ~0  // all bits set
#define DROPSDebugC 0 

/// \brief Throws an error upon a failed assertion.
///
/// \param a The assertion; must be convertible to bool.
/// \param b The object to be thrown, if a==false.
/// \param c The debugging-class, to which this assertion belongs.
/// \remarks Using a macro ensures, that the compiler (and optimizer)
/// is not confused by the template-functions and classes used for error
/// reporting in portions of the code that shall not be debugged.
#if DROPSDebugC
#  define Assert(a,b,c) (_Assert((a),(b),(c)))
#else
#  define Assert(a,b,c) ((void)0)
#endif

/// \brief Conditionally write a message to the debugging stream.
///
/// The condition will be checked if any debugging-class is active. If the
/// condition is true, the message will be written to the debug-stream.
/// \param a The message to be written. This must be an expression suitable
///     for writing after cdebug <<.
/// \param b The condition, under which the message is written; must be convertible
///     to bool.
#if DROPSDebugC
#  define Comment(a,b) do { if ((b) & DROPSDebugC) cdebug << a; } while (false)
#else
#  define Comment(a,b) ((void)0)
#endif

/// Shut up gcc to not warn about certain unused function-parameters.
#ifdef __GNUC__
#  define __UNUSED__ __attribute__((__unused__))
#else
#  define __UNUSED__
#endif

/// \name Macros for valarray-derivatives.
/// Several of the numerical classes, e.g. VectorCL, QuadbaseCL
/// LocalP2CL, etc, are derived from valarray. To take advantage of expression
/// template mechanisms some constructors and copy-assignment ops must be
/// defined. This repititive task is simplified by the following makros.
//@{
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
    Assert( this->size()==thebase_type( v).size(),                                       \
            #theClass #theOp ": incompatible dimensions", DebugNumericC);                \
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
/// \brief Call this macro in the definition of a class that is derived
///     from valarray.
///
/// \param theClass Name of the derived class.
/// \param theT Type of the components of the valarray.
/// \param thebase_type Name of the immidiate base-class.
#define DROPS_DEFINE_VALARRAY_DERIVATIVE(theClass, theT, thebase_type)     \
/*The expression template constructor*/                                    \
DROPS_EXP_TEMPL_CONSTR_FOR_VALARRAY_DERIVATIVE(theClass, thebase_type)     \
/*assignment and computed assignment*/                                     \
DROPS_ASSIGNMENT_OPS_FOR_VALARRAY_DERIVATIVE(theClass, theT, thebase_type)
//@}


/// \brief Check, if a value is in a sequence.
///
/// Returns true, iff value is in [beg, end).
template <class In, class T>
inline bool is_in( In beg, In end, const T& value)
{
    return std::find(beg,end,value) != end;
}


/// \brief Check, if a predicate holds anywhere in a sequence.
///
/// Returns true, iff there is v in [beg, end) with p( v) == true.
template <class In, class Pred>
inline bool is_in_if( In beg, In end, Pred p )
{
    return std::find_if(beg,end,p) != end;
}


/// \brief Functor, that converts a reference to a pointer.
///
/// Useful for some STL-like algorithms.
template <class T>
class ref_to_ptr : public std::unary_function<T&, T*>
{
  public:
    T* operator() (T& arg) const { return static_cast<T*>(&arg); }
};


/// \brief Base class for all classes that DROPS throws as exceptions.
///
/// Classes should derive their own (hopefully more powerful) error-class
/// and donate meaningful error-messages. The default error handler only
/// prints the error message and abort.
class DROPSErrCL
{
  protected:
    std::string _ErrMesg;

  public:
    DROPSErrCL() : _ErrMesg("") {}
    DROPSErrCL(const std::string& mesg) : _ErrMesg(mesg) {}
    virtual ~DROPSErrCL() {}

    /// Override this to inform the user about details of the error-condition.
    virtual std::ostream& what  (std::ostream&) const;
    /// Lets you provide your own error-handler.
    virtual void handle() const;
};


/// Used by the Assert macro.
//@{
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
//@}


/// \brief Provides a unique identifier for an object.
///
/// We use the template argument to specify the class, whose objects will
/// carry an Id. Users of the code should never need to construct these
/// objects themselves.
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

    /// Used by MakeConsistentNumbering().
    static void ResetCounter(Ulint i= 0) { _Counter= i; }

    bool operator == (const IdCL<type>& Id) const
        { return Id._Identity == _Identity; }
    bool operator != (const IdCL<type>& Id) const { return !(*this==Id); }
    bool operator <  (const IdCL<type>& Id) const
        { return _Identity < Id._Identity; }
};

template <class type> Ulint IdCL<type>::_Counter = 0;


/// Get to know how fast DROPS is !  :-)
class TimerCL
{
  private:
    rusage _t_begin, _t_end;
    double _time;

  public:
    TimerCL(double time=0.) { Reset(time); }

    /// Start, stop, and read timer
    //@{
    double GetTime() const     { return _time; }
    void Reset(double time= 0) { _time= time; Start(); }
    void Start()               { getrusage(RUSAGE_SELF,&_t_begin); }
    void Stop()
    {
        getrusage(RUSAGE_SELF,&_t_end);
        _time+= (_t_end.ru_utime.tv_sec - _t_begin.ru_utime.tv_sec)
              + (_t_end.ru_stime.tv_sec - _t_begin.ru_stime.tv_sec)
              + double(_t_end.ru_utime.tv_usec - _t_begin.ru_utime.tv_usec
                      +_t_end.ru_stime.tv_usec - _t_begin.ru_stime.tv_usec)/1000000; }
    //@}
};


/// \brief Functor to select the second component of a std::pair-like type.
///
/// This is needed, as the C++-standard comittee deemed selectors for
/// pairs unneccessary.
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


/// \brief Iterator for a sequence of objects that is given as a sequence of pointers
///     to these objects.
///
/// This is a random access iterator. The iterator_traits of the standard-library work
/// with this class as the neccessary typedefs are defined.
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


/// \brief Deal with const-qualification in template-metaprogramming.
///
/// For a given type T, stripped_type is T with a possible outermost const removed,
/// const_type adds a const-qualifier if T did not have one.
//@{
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
//@}

/// \brief Deal with value_type of containers in template-metaprogramming.
///
/// For a given type container-type T, value_type is the type of the values stored in T
/// This works for pointers and arrays, too.
//@{
template<class T>
  struct ValueHelperCL
{
    typedef typename T::value_type value_type;
};

template<class T>
  struct ValueHelperCL<T*>
{
    typedef T value_type;
};
template<class T, size_t S>
  struct ValueHelperCL<T[S]>
{
    typedef T value_type;
};
//@}

} // end of namespace DROPS

#endif
