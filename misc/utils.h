/// \file
/// \brief Useful stuff that fits nowhere else.

#ifndef DROPS_UTILS_H
#define DROPS_UTILS_H


// #include <limits> ///< \todo Do we have limits with gcc-snapshots or SGI-CC?
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#if defined(DROPS_WIN)
#  include <windows.h>
#  undef max
#  undef min
#else
#  include <sys/time.h>
#  include <sys/resource.h>
#  include <sys/stat.h>
#  include <sys/types.h>
#endif
#include <cmath>

#ifndef M_PI
# ifdef DROPS_WIN
#  define _USE_MATH_DEFINES
# endif
# include <math.h>
#endif

#ifdef _OPENMP
#  include <omp.h>  // for timing
#endif

#ifdef DROPS_WIN
double cbrt(double arg);
#endif

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

/// Master process
#ifdef _PAR
#  define Drops_MasterC 0
#  define MASTER (DROPS::ProcCL::IamMaster())
#  define IF_MASTER if (MASTER)
#  define IF_NOT_MASTER if (!MASTER)
/// Uncomment the following line to use C++-interface of MPI
//#   define _MPICXX_INTERFACE
#else
#  define MASTER true
#  define IF_MASTER
#  define IF_NOT_MASTER if (false)
#endif

/// \name Code-groups for debugging.
/// \brief Constants that group the code for debugging and error-reporting.
//@{
#define DebugContainerC      1
#define DebugRefineEasyC     2
#define DebugRefineHardC     4
#define DebugNumericC        8
#define DebugUnknownsC      16
#define DebugNoReuseSparseC 32
#define DebugParallelC      64
#define DebugParallelHardC 128
#define DebugParallelNumC  256
#define DebugLoadBalC      512
#define DebugDiscretizeC  1024
#define DebugSubscribeC   2048
#define DebugOutPutC      4096
//@}

/// The stream for dedug output.
/// In parallel mode, the proc number is printed in front of the message
#ifndef _PAR
#  define cdebug std::cout
#else
#  define cdebug std::cout << "["<<ProcCL::MyRank()<<"]: "
#endif

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
/// In parallel mode only the master-process will write this comment.
/// If the comment should be written by any process use AllComment.
/// \param a The message to be written. This must be an expression suitable
///     for writing after cdebug <<.
/// \param b The condition, under which the message is written; must be convertible
///     to bool.
#if DROPSDebugC
#  define AllComment(a,b) do { if ((b) & DROPSDebugC) cdebug << a; } while (false)
#  define Comment(a,b) do { if ((b) & DROPSDebugC) IF_MASTER cdebug << a; } while (false)
#else
#  define AllComment(a,b) ((void)0)
#  define Comment(a,b) ((void)0)
#endif

/// Shut up gcc to not warn about certain unused function-parameters.
#ifdef __GNUC__
#  define __UNUSED__ __attribute__((__unused__))
#else
#  define __UNUSED__
#endif

/// \brief Select how to handle negative norms
#define DROPS_ABORT_ON_NEG_SQ_NORM

/// \brief Constant for zero squared norm
const double ZeroNormC = 1.0e-32;

/// \brief Makro to handle negative squared norm
///
/// It may happen, that a squared norm is negative due to rounding errors while
/// computing the norm of a distributed vector. This makro set the norm to
/// ZeroNormC or throws an exception.
#ifdef DROPS_ABORT_ON_NEG_SQ_NORM
#define DROPS_Check_Norm(a,b) if ((a)<0.) { cdebug << "Norm is "<<(a)<< std::endl; throw DROPSErrCL((b)); }
#else
#define DROPS_Check_Norm(a,b) if ((a)<0.) (a)= ZeroNormC
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


/// \brief Get the address of the first element in a valarray
///
/// ("&x[0]" doesn't work, because "operator[] const" only returns a value)
/// \todo (merge)  Addr() functions in 'misc/utils.h'?
template <typename T>
  inline const T*
  Addr(const std::valarray<T>& x)
{
    return &(const_cast<std::valarray<T>&>(x)[0]);
}

/// \brief Get the address of the first element in a valarray
template <typename T>
  inline T*
  Addr(std::valarray<T>& x)
{
    return &(x[0]);
}

template <typename T>
  inline const T*
  Addr(const std::vector<T>& x)
{
    return &(const_cast<std::vector<T>&>(x)[0]);
}

/// \brief Get the address of the first element in a vector
template <typename T>
  inline T*
  Addr(std::vector<T>& x)
{
    return &(x[0]);
}

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


/// \brief Iterate through a STL-container and do an operation if a condition holds
template <class In, class Op, class Pred>
inline void for_each_if( In beg, In end, Op f, Pred p )
{
    while (beg!=end) { if (p(*beg)) f(*beg); ++beg; }
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
        if (!assertion)
            throw exc;
}


template <class A>
inline void
_Assert(A assertion, const char* msg, Uint DebugLevel=~0)
{
    if (DebugLevel&DROPSDebugC)
        if (!assertion)
            throw DROPSErrCL(msg);
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


/// \brief Get to know how fast DROPS is !  :-)
/// Time measurement is done by getrusage if OpenMP is not enabled, otherwise
/// use OpenMP to determine time.
class TimerCL
{
  private:
#ifdef _OPENMP
    double _t_begin;
#else
    rusage _t_begin, _t_end;
#endif
    double _time;


  public:
    TimerCL(double time=0.) { Reset(time); }

    /// Start, stop, and read timer
    //@{
    double GetTime() const     { return _time; }
    void Reset(double time= 0) { _time= time; Start(); }
    void Start()               {
#ifndef _OPENMP
        getrusage(RUSAGE_SELF,&_t_begin);
#else
        _t_begin= omp_get_wtime();
#endif
        }
    void Stop()
    {
#ifndef _OPENMP
        getrusage(RUSAGE_SELF,&_t_end);
        _time+= (_t_end.ru_utime.tv_sec - _t_begin.ru_utime.tv_sec)
              + (_t_end.ru_stime.tv_sec - _t_begin.ru_stime.tv_sec)
              + double(_t_end.ru_utime.tv_usec - _t_begin.ru_utime.tv_usec
                      +_t_end.ru_stime.tv_usec - _t_begin.ru_stime.tv_usec)/1000000;
#else
        _time+= omp_get_wtime()-_t_begin;
#endif
        }
    //@}
};


/// \brief Represents the permutation i-->p[i] on [0, ..., p.size()).
typedef std::vector<size_t> PermutationT;

/// \brief Compute the inverse permutation of p, id est pi[p[i]] == i for all i.
PermutationT
invert_permutation (const PermutationT& p);


/// \brief Return the first iterator i from [begin, end) where l(*j, *i) == false for all j in [begin, end) excluding i.
/// For an empty sequence end is returned.
/// \todo (merge) What is the difference to std::min_element?
template <class Iterator, class Cmp>
  Iterator
  arg_min (Iterator begin, Iterator end, Cmp l= std::less<typename std::iterator_traits<Iterator>::value_type>())
{
    if (begin == end) return end;

    typename std::iterator_traits<Iterator>::value_type m= *begin;
    Iterator am= begin;
    while (++begin != end)
        if ( l( *begin, m)) {
            m= *begin;
            am= begin;
        }
    return am;
}


/// \brief Return the first iterator i from [begin, end) where *i<=*j for all j in [begin, end).
/// For an empty sequence end is returned.
/// \todo (merge) What is the difference to std::min_element?
template <class Iterator>
  Iterator
  arg_min (Iterator begin, Iterator end)
{
    return arg_min( begin, end,
        std::less<typename std::iterator_traits<Iterator>::value_type>());
}

/// \brief Output [begin, end) to out, separated by newline.
template <class Iterator>
void
inline seq_out (Iterator begin, Iterator end, std::ostream& out)
{
    for (; begin != end; ++begin) out << *begin << '\n';
}

/// \brief Output obj via operator<<  to a file filename.
///
/// The filename and an optional name are reported on std::cout.
template <class StreamableT>
void
WriteToFile (const StreamableT& obj, std::string filename , std::string name= std::string())
{
    std::ofstream mystream( filename.c_str());
    mystream.precision( 18);
    if (!mystream) {
        std::cout << filename << std::endl;
        throw DROPSErrCL( "WriteToFile: error while opening file\n");
    }
    std::cout << "Writing to file \"" << filename << "\".    Description: " << name << '\n';
    mystream << obj << std::flush;
    if (!mystream)
        throw DROPSErrCL( "WriteToFile: write failed\n");
}

/// \brief Functor to select the second component of a std::pair-like type.
///
/// This is needed, as the C++-standard committee deemed selectors for
/// pairs unnecessary.
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

/// \brief Predicate that compares a std::pair-like type by its first
///     component only.
template <class Pair>
struct less1st: public std::binary_function<Pair, Pair, bool>
{
  bool operator() (const Pair& x, const Pair& y) const {
    return x.first < y.first;
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

/// \brief Create a directory
int CreateDirectory(std::string path);

/// \brief Remove a file
int DeleteFile(std::string file);

/// \brief stream buffer without output on screen (as writing to /dev/null)
class NullStreambufCL: public std::streambuf {
  public:
	NullStreambufCL() {}

  protected:
	/// \brief Usually this streambuf member is used to write output to some physical device. Here the output is just ignored.
    virtual int_type overflow( int_type c)
    { return c; }
};

/// \brief Mute and restore standard output streams
class MuteStdOstreamCL
{
private:
    std::streambuf *bout_, *berr_, *blog_;
    NullStreambufCL devnull_;

public:
    MuteStdOstreamCL()
    : bout_(std::cout.rdbuf()), berr_(std::cerr.rdbuf()), blog_(std::clog.rdbuf()) {}
    /// Mute given stream
    void Mute( std::ostream& os) { os.rdbuf( &devnull_); }
    /// Mute std::cout, std::cout, std::clog
    void Mute() { Mute(std::cout); Mute(std::clog); }
    /// Recover behavior of std::cout, std::cout, std::clog prior construction of this object
    void Recover() const { std::cout.rdbuf(bout_); std::cerr.rdbuf(berr_); std::clog.rdbuf(blog_); }
};

} // end of namespace DROPS

#ifdef _PAR
#  ifndef _MPICXX_INTERFACE
#    define MPICH_SKIP_MPICXX
#  endif
#  pragma GCC system_header  // Suppress warnings from mpi.h
#  include <mpi.h>
#endif

#endif
