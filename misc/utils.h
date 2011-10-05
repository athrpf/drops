/// \file utils.h
/// \brief Useful stuff that fits nowhere else.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Martin Horsky, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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
#include <map>
#if defined(DROPS_WIN)
#  include <windows.h>
#  undef max
#  undef min
#  include <unordered_map>
#  define DROPS_STD_UNORDERED_MAP std::unordered_map
#else
#  include <sys/time.h>
#  include <sys/resource.h>
#  include <sys/stat.h>
#  include <sys/types.h>
#endif
#include <cmath>
#if __GNUC__ >= 4
#    include <tr1/unordered_map>
#    define DROPS_STD_UNORDERED_MAP std::tr1::unordered_map
#endif

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

#include <cstddef>

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
#ifndef DROPSDebugC 
  #define DROPSDebugC 0
#endif  

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

/// \brief Predicate, that compares a std::pair-like type by its first
///     component only.
template <class Pair>
struct less1st: public std::binary_function<Pair, Pair, bool>
{
  bool operator() (const Pair& x, const Pair& y) const {
    return x.first < y.first;
  }
};

/// \brief Predicate, that compares a pointers by the values, to which they point.
template <class PtrT>
struct less_by_ptr: public std::binary_function<PtrT, PtrT, bool>
{
  bool operator() (PtrT x, PtrT y) const {
    return *x < *y;
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


/// \brief Provide begin()-iterator uniformly for standard-containers, std::valarray and arrays.
///@{
template <class ContainerT>
class SequenceTraitCL
{
  public:
    typedef ContainerT container_type;

    typedef typename ContainerT::iterator             iterator;
    typedef typename ContainerT::const_iterator const_iterator;

    static       iterator       begin (container_type& c)       { return c.begin(); }
    static const_iterator const_begin (const container_type& c) { return c.begin(); }
};

template <class T>
class SequenceTraitCL<std::valarray<T> >
{
  public:
    typedef std::valarray<T> container_type;

    typedef       T*       iterator;
    typedef const T* const_iterator;

    static       iterator       begin (container_type& c)       { return Addr( c); }
    static const_iterator const_begin (const container_type& c) { return Addr( c); }
};

template <class T>
  class GridFunctionCL; ///< forward declaration

template <class T>
class SequenceTraitCL<GridFunctionCL<T> >
{
  public:
    typedef GridFunctionCL<T> container_type;

    typedef       T*       iterator;
    typedef const T* const_iterator;

    static       iterator       begin (container_type& c)       { return Addr( c); }
    static const_iterator const_begin (const container_type& c) { return Addr( c); }
};

template <class T, size_t S>
class SequenceTraitCL<T[S]>
{
  public:
    typedef       T*       iterator;
    typedef const T* const_iterator;

    static       iterator       begin (T c[S])       { return c; }
    static const_iterator const_begin (const T c[S]) { return c; }
};

///\brief returns the begin of c for container-types, valarrays and arrays.
template <class ContainerT>
  inline typename SequenceTraitCL<ContainerT>::iterator
  sequence_begin (ContainerT& c) { return SequenceTraitCL<ContainerT>::begin( c); }

///\brief returns the begin of c for const-container-types, const-valarrays and const-arrays.
template <class ContainerT>
  inline typename SequenceTraitCL<ContainerT>::const_iterator
  sequence_begin (const ContainerT& c) { return SequenceTraitCL<ContainerT>::const_begin( c); }
///@}

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

/// \brief Reversal of the byte order (change from little to big endian decoding)
void reverseByteOrder(int size,char field[]);


//@{ used by error marker
class TetraCL;
typedef std::pair<const TetraCL*, double> Err_PairT;
typedef std::vector<Err_PairT> Err_ContCL;

struct AccErrCL :public std::binary_function<double, const Err_PairT, double>
{
    double operator() (double init, const Err_PairT& ep) const
        { return init + ep.second;}
};

struct Err_Pair_GTCL :public std::binary_function<const Err_PairT, const Err_PairT, bool>
{
    bool operator() (const Err_PairT& ep0, const Err_PairT& ep1)
    { return ep0.second > ep1.second; }
};
//@}

/// \brief Transforms a std::map into a vector of pairs
/** This can be used to parallelize a loop over all
    elements in a map.
*/
template <typename T1, typename T2>
std::vector<std::pair<T1,T2> > Map2Vec( const std::map<T1,T2>& Map)
{
    std::vector<std::pair<T1,T2> > vec(Map.size());
    size_t pos=0;
    for ( typename std::map<T1,T2>::const_iterator it=Map.begin(); it!=Map.end(); ++it, pos++)
        vec[pos]= std::pair<T1,T2>(it->first, it->second);
    return vec;
}

#if __GNUC__ >= 4 && !defined(__INTEL_COMPILER)
template <typename T1, typename T2>
std::vector<std::pair<T1,T2> > Map2Vec( const std::tr1::unordered_map<T1,T2>& Map)
{
    std::vector<std::pair<T1,T2> > vec(Map.size());
    size_t pos=0;
    for ( typename std::tr1::unordered_map<T1,T2>::const_iterator it=Map.begin(); it!=Map.end(); ++it, pos++)
        vec[pos]= std::pair<T1,T2>(it->first, it->second);
    return vec;
}
#endif

inline bool
logical_xor (bool a, bool b)
{
    return ( a || b) && !(a && b);
}

///\brief The sign of the argument in \f$\{-1,0,+1\}\f$
inline byte sign (double d)
{
    return d > 0. ? 1 : (d < 0. ? -1 : 0);
}

} // end of namespace DROPS


#ifdef DROPS_WIN
///\brief Assignement of slice array is missing in VS Compi
template<typename T>
inline std::slice_array<T>&
std::slice_array<T>::operator=(const slice_array<T>& a)
{
    for (size_t i = 0; i < a.size(); ++i){
        this->_Myptr[i * this->stride()] = a._Myptr[i *a.stride()];
	}
	return *this;
}
#endif

#ifdef _PAR
#  ifndef _MPICXX_INTERFACE
#    define MPICH_SKIP_MPICXX
#  endif
#  ifndef DROPS_WIN
#    pragma GCC system_header  // Suppress warnings from mpi.h
#  endif
#  include <mpi.h>
#  ifdef _HYPRE
#    include <HYPRE.h>
#    include <HYPRE_IJ_mv.h>
#    include <HYPRE_parcsr_ls.h>
#  endif
#endif

#ifndef DROPS_WIN
#  pragma GCC system_header // Suppress warnings from boost
#endif
# include <boost/property_tree/ptree.hpp>
# include <boost/property_tree/exceptions.hpp>
# include <boost/property_tree/json_parser.hpp>


#endif
