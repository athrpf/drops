//**************************************************************************
// File:     container.h                                                   *
// Content:  simple array in STL style                                     *
//           level-based array                                             *
// Author:   Joerg Peters, Volker Reichelt, IGPM RWTH Aachen               *
// Version:  0.2                                                           *
// History:                                                                *
//                                                                         *
// Remarks:                                                                *
//                                                                         *
//**************************************************************************


#ifndef DROPS_CONTAINER_H
#define DROPS_CONTAINER_H

#include <vector>
#include <list>
#include <cmath>
#include <iostream>
#include "misc/utils.h"

namespace DROPS
{

//**************************************************************************
// Class:    DMatrixCL                                                     *
// Purpose:  dynamical storage for 2D data                                 *
// Remarks:  just the bare minimum                                         *
//**************************************************************************

template <typename T>
class DMatrixCL
{
  private:
    size_t Cols_;
    T* Array_;

  public:
    DMatrixCL(size_t row, size_t col) : Cols_(col), Array_(new T[row*col]) {}
    ~DMatrixCL() { delete[] Array_; }

    T& operator() (size_t row, size_t col)       { return Array_[row*Cols_+col]; }
    T  operator() (size_t row, size_t col) const { return Array_[row*Cols_+col]; }
};


template <class T, Uint _Size>
  class SArrayCL;
template <class T, Uint _Size>
  class SBufferCL;
template <Uint _Size>
  class SVectorCL;
template <Uint _Rows, Uint _Cols>
  class SMatrixCL;
template <class T, Uint _Size>
  inline bool
  operator==(const SArrayCL<T, _Size>&, const SArrayCL<T, _Size>&);
template <class T, Uint _Size>
  inline bool
  operator==(const SBufferCL<T, _Size>&, const SBufferCL<T, _Size>&);

typedef SVectorCL<2> Point2DCL;
typedef SVectorCL<3> Point3DCL;
typedef SVectorCL<4> BaryCoordCL;


//**************************************************************************
// Class:    SArrayCL                                                      *
// Purpose:  an array that remembers its size                              *
// Remarks:  All functions are inline, should be as fast as a "bare" array *
//**************************************************************************

template <class T, Uint _Size>
class SArrayCL
{
private:
    T Array[_Size];

public:
    typedef       T*       iterator;
    typedef const T* const_iterator;
    typedef       T&       reference;
    typedef const T& const_reference;
    typedef       T  value_type;

//    SArrayCL() {}
    explicit           SArrayCL(T val= T())       { std::fill_n(Array+0, _Size, val); }
    template<class In> SArrayCL(In start)         { std::copy(start, start+_Size, Array+0); }
    template<class In> SArrayCL(In start, In end) { std::copy(start, end, Array+0); }
    // Default copy-ctor, assignment operator, dtor

    iterator       begin     ()             { return static_cast<T*>(Array); }
    const_iterator begin     ()       const { return static_cast<const T*>(Array); }
    iterator       end       ()             { return Array+_Size; }
    const_iterator end       ()       const { return Array+_Size; }
    reference      operator[](Uint i)       { Assert(i<_Size, DROPSErrCL("SArrayCL::operator[]: wrong index"), DebugContainerC);
                                              return Array[i]; }
    value_type     operator[](Uint i) const { Assert(i<_Size, DROPSErrCL("SArrayCL::operator[]: wrong index"), DebugContainerC);
                                              return Array[i]; }
    Uint           size      ()       const { return _Size; }

    friend bool operator==<>(const SArrayCL&, const SArrayCL&);
};

template <class T, Uint _Size>
  inline bool
  operator==(const SArrayCL<T, _Size>& a0, const SArrayCL<T, _Size>& a1)
{
    for (Uint i=0; i<_Size; ++i)
        if (a0[i] != a1[i]) return false;
    return true;
}

//**************************************************************************
// Class:    SBufferCL                                                     *
// Purpose:  A buffer or ring-buffer of fixed size with wrap-around indices*
//**************************************************************************
template <class T, Uint _Size>
class SBufferCL
{
private:
    T Array[_Size];
    int Front;

public:
    typedef       T&       reference;
    typedef const T& const_reference;
    typedef       T  value_type;

//    SBufferCL() {}
    explicit           SBufferCL(T val= T())      { std::fill_n(Array+0, _Size, val); Front= 0; }
    template<class In> SBufferCL(In start)         { std::copy(start, start+_Size, Array+0); Front= 0; }
    template<class In> SBufferCL(In start, In end) { std::copy(start, end, Array+0); Front= 0; }
    SBufferCL(const SBufferCL& b) {
        std::copy( b.Array+0, b.Array+_Size, Array+0);
        Front= b.Front; }
    SBufferCL& operator=(const SBufferCL& b) {
        if (&b!=this) {
            std::copy( b.Array+0, b.Array+_Size, Array+0);
            Front= b.Front; }
        return *this; }
    // Default dtor

    reference  operator[](int i) {
        if (i<0) i+= _Size;
        else if (i >= static_cast<int>( _Size)) i-=_Size;
        Assert( (i>=0 && i<static_cast<int>( _Size)), DROPSErrCL("SBufferCL::operator[]: wrong index"), DebugContainerC);
        return Array[(i+Front < static_cast<int>( _Size)) ? i+Front : i+Front-_Size]; }
    value_type operator[](int i) const {
        if (i<0) i+= _Size;
        else if (i >= static_cast<int>( _Size)) i-=_Size;
        Assert( (i>=0 && i<_Size), DROPSErrCL("SBufferCL::operator[]: wrong index"), DebugContainerC);
        return Array[(i+Front < static_cast<int>( _Size)) ? i+Front : i+Front-_Size]; }
    Uint size() const { return _Size; }

    // first entry is overwritten.
    void push_back(const T& t) {
        Array[Front]= t;
        if (++Front == static_cast<int>( _Size)) Front-= _Size; }
    // for positive i, front element moves to the end.
    void rotate (int i= 1) {
        Assert(i<static_cast<int>( _Size), DROPSErrCL("SBufferCL::rotate: wrong index"), DebugContainerC);
        Front+= i;
        if (Front >= static_cast<int>( _Size)) Front-= _Size;
        else if (Front < 0) Front+= _Size; }
        

    friend bool operator==<>(const SBufferCL&, const SBufferCL&);
};

template <class T, Uint _Size>
  inline bool
  operator==(const SBufferCL<T, _Size>& a0, const SBufferCL<T, _Size>& a1)
{
    if (a0.Front != a1.Front) return false;
    for (Uint i=0; i<_Size; ++i)
        if (a0.Array[i] != a1.Array[i]) return false;
    return true;
}

//**************************************************************************
// Class:    SVectorCL                                                     *
// Purpose:  primitive vector class with templates for short vectors       *
// Remarks:  Many optimizations are possible.                              *
//**************************************************************************

template <Uint _Size>
class SVectorCL : public SArrayCL<double,_Size>
{
public:
    SVectorCL() {}
    explicit           SVectorCL(double val)       : SArrayCL<double,_Size>(val)       {}
//    SVectorCL(const SliceArrayCL& s)
//        { for (Uint i=0; i!=_Size; ++i) (*this)[i]= s[i]; }
    template<class In> SVectorCL(In start)         : SArrayCL<double,_Size>(start)     {}
    template<class In> SVectorCL(In start, In end) : SArrayCL<double,_Size>(start,end) {}

//    SVectorCL& operator= (const SliceArrayCL&);
    SVectorCL& operator+=(const SVectorCL&);
    SVectorCL& operator-=(const SVectorCL&);
    SVectorCL& operator*=(double);
    SVectorCL& operator/=(double);
    
    // komponentenweise Operatoren
    SVectorCL& operator*=(const SVectorCL&);
    SVectorCL& operator/=(const SVectorCL&);

    double norm_sq() const;
    double norm()    const { return ::sqrt(norm_sq()); }
};

/*
template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator=(const SliceArrayCL& s)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]= s[i];
    return *this;
}
*/

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator+=(const SVectorCL& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]+= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator-=(const SVectorCL<_Size>& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]-= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator*=(const SVectorCL<_Size>& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]*= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator/=(const SVectorCL<_Size>& v)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]/= v[i];
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator*=(double s)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]*= s;
    return *this;
}

template <Uint _Size>
SVectorCL<_Size>&
SVectorCL<_Size>::operator/=(double s)
{
    for (Uint i=0; i!=_Size; ++i) (*this)[i]/= s;
    return *this;
}

template <Uint _Size>
double SVectorCL<_Size>::norm_sq() const
{
    double x = 0.0;
    for (Uint i=0; i<_Size; ++i) x += (*this)[i]*(*this)[i];
    return x;
}

template <Uint _Size>
SVectorCL<_Size> BaryCenter(const SVectorCL<_Size>& v1, const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i] = .5 * (v1[i] + v2[i]);
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> ConvexComb (double a,
                             const SVectorCL<_Size>& v1,
                             const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i] = (1.0-a)*v1[i] + a*v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator+(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] + v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator-(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] - v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator-(const SVectorCL<_Size>& v1)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= -v1[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator*(double d, const SVectorCL<_Size>& v)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= d * v[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator*(const SVectorCL<_Size>& v, double d)
{
    return d*v;
}

template <Uint _Size>
double inner_prod(const SVectorCL<_Size>& v1, const SVectorCL<_Size>& v2)
{
    double ret= 0.0;
    for (Uint i=0; i <_Size; ++i) ret+= v1[i]*v2[i];
    return ret;
}

template <Uint _Size>
SVectorCL<_Size> operator/(const SVectorCL<_Size>& v, double d)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= v[i]/d;
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator*(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] * v2[i];
    return tempv;
}

template <Uint _Size>
SVectorCL<_Size> operator/(const SVectorCL<_Size>& v1,
                           const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= v1[i] / v2[i];
    return tempv;
}

template <Uint _Size>
bool operator<(const SVectorCL<_Size>& v1,
               const SVectorCL<_Size>& v2)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) if(!( v1[i] < v2[i]) ) return false;
    return true;
}

using ::sqrt;

template <Uint _Size>
SVectorCL<_Size> sqrt(const SVectorCL<_Size>& v)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= sqrt(v[i]);
    return tempv;
}

using ::fabs;

template <Uint _Size>
SVectorCL<_Size> fabs(const SVectorCL<_Size>& v)
{
    SVectorCL<_Size> tempv;
    for (Uint i=0; i<_Size; ++i) tempv[i]= fabs(v[i]);
    return tempv;
}


template <Uint _Size>
std::ostream& operator<<(std::ostream& os, const SVectorCL<_Size>& v)
{
//    os << v.size() << "    ";
    for (Uint i=0; i<v.size(); ++i)
        os << v[i] << ' ';
    return os;
}


inline void
cross_product(Point3DCL& res, const Point3DCL& v0, const Point3DCL& v1)
// res= v0 x v1
{
    res[0]= v0[1]*v1[2] - v0[2]*v1[1];
    res[1]= v0[2]*v1[0] - v0[0]*v1[2];
    res[2]= v0[0]*v1[1] - v0[1]*v1[0];
}

// std_basis<n>(0)==Null, std_basis<n>(i)[j]==Delta_i-1_j
template <Uint _Size>
inline SVectorCL<_Size> std_basis(Uint i)
{
    SVectorCL<_Size> ret(0.);
    if (i>0) ret[i-1]= 1.;
    return ret;
}


template <Uint _Rows, Uint _Cols>
class SMatrixCL : public SVectorCL<_Rows*_Cols>
{
  public:
    typedef SVectorCL<_Rows*_Cols> _vec_base;

    SMatrixCL() {}
    explicit           SMatrixCL(double val)       : SVectorCL<_Rows*_Cols>(val)       {}
    template<class In> SMatrixCL(In start)         : SVectorCL<_Rows*_Cols>(start)     {}
    template<class In> SMatrixCL(In start, In end) : SVectorCL<_Rows*_Cols>(start,end) {}

// Schreib- & Lesezugriff
    double& operator() (int row, int col)       { return (*this)[row*_Cols+col]; }// Matrix(i,j)
    double  operator() (int row, int col) const { return (*this)[row*_Cols+col]; }
//    SliceArrayCL operator() (const SliceCL& sl)
//        { return SliceArrayCL(this->begin(), sl); }

// Zuweisung & Co.
    SMatrixCL& operator+=(const SMatrixCL&);                // Matrix=Matrix+Matrix'
    SMatrixCL& operator-=(const SMatrixCL&);                // Matrix=Matrix-Matrix'
    SMatrixCL& operator*=(double s);                               // Matrix = c * Matrix
    SMatrixCL& operator/=(double s);                               // Matrix = Matrix'/c

// Dimensionen feststellen
    Uint num_rows() const { return _Rows; }                        // Zeilenzahl
    Uint num_cols() const { return _Cols; }                        // Spaltenzahl
};

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator+=(const SMatrixCL<_Rows, _Cols>& m)
{
    *static_cast<_vec_base*>(this)+= *static_cast<const _vec_base*>(&m);
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator-=(const SMatrixCL<_Rows, _Cols>& m)
{
    *static_cast<_vec_base*>(this)-= *static_cast<const _vec_base*>(&m);
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator*=(double d)
{
    *static_cast<_vec_base*>(this)*= d;
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>&
SMatrixCL<_Rows, _Cols>::operator/=(double d)
{
    *static_cast<_vec_base*>(this)/= d;
    return *this;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator+(const SMatrixCL<_Rows, _Cols>& m1, const SMatrixCL<_Rows, _Cols>& m2)
{
    SMatrixCL<_Rows, _Cols> ret;
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m1)
         +*static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m2);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator-(const SMatrixCL<_Rows, _Cols>& m1, const SMatrixCL<_Rows, _Cols>& m2)
{
    SMatrixCL<_Rows, _Cols> ret;
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m1)
         -*static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m2);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator-(const SMatrixCL<_Rows, _Cols>& m)
{
    SMatrixCL<_Rows, _Cols> ret;
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = -*static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator*(double d, const SMatrixCL<_Rows, _Cols>& m)
{
    SMatrixCL<_Rows, _Cols> ret;
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = d**static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator*(const SMatrixCL<_Rows, _Cols>& m, double d)
{
    SMatrixCL<_Rows, _Cols> ret;
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m)*d;
    return ret;
}

template<Uint _Rows, Uint _Cols>
SMatrixCL<_Rows, _Cols>
operator/(const SMatrixCL<_Rows, _Cols>& m, double d)
{
    SMatrixCL<_Rows, _Cols> ret;
    *static_cast<typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&ret)
        = *static_cast<const typename SMatrixCL<_Rows, _Cols>::_vec_base*>(&m)/d;
    return ret;
}

template<Uint _RowsL, Uint _ColsR, Uint _Dim>
SMatrixCL<_RowsL, _ColsR>
operator*(const SMatrixCL<_RowsL, _Dim>& m1, const SMatrixCL<_Dim, _ColsR>& m2)
{
    SMatrixCL<_RowsL, _ColsR> ret(0.0);
    for (Uint row=0; row!=_RowsL; ++row)
        for (Uint col=0; col!=_ColsR; ++col)
            for (Uint i=0; i!=_Dim; ++i)
                ret(row, col)+= m1(row, i)*m2(i, col);
    return ret;
}

template<Uint _Rows, Uint _Cols>
SVectorCL<_Rows>
operator*(const SMatrixCL<_Rows, _Cols>& m, const SVectorCL<_Cols>& v)
{
    SVectorCL<_Rows> ret(0.0);
    for (Uint row=0; row!=_Rows; ++row)
        for (Uint i=0; i!=_Cols; ++i)
                ret[row]+= m(row, i)*v[i];
    return ret;
}

template <Uint _Rows, Uint _Cols>
std::ostream& operator << (std::ostream& os, const SMatrixCL<_Rows, _Cols>& m)
{
    const Uint M = m.num_rows();
    const Uint N = m.num_cols();

    os << M << ' ' << N << '\n' ;

    for (Uint row=0; row<M; ++row)
    {
        os << "  ";
        for (Uint col=0; col<N; ++col)
            os << m(row, col) << ' ';
        os << '\n';
    }
    return os;
}


//**************************************************************************
// Class:    GlobalListCL                                                  *
// Purpose:  List of lists per level                                       *
// Remarks:  efficient random access for levels                            *
//**************************************************************************

template <class T>
class GlobalListCL
{
  public:
    typedef std::list<T>                                  LevelCont;
    typedef typename LevelCont::iterator                  LevelIterator;
    typedef typename LevelCont::const_iterator            const_LevelIterator;
    // TODO: If needed, make this an efficient random access iterator using _LevelStarts
    typedef typename std::list<LevelCont>::iterator       iterator;
    typedef typename std::list<LevelCont>::const_iterator const_iterator;

  private:
    Uint                 _LevelBegin;
    std::list<LevelCont> _Data;
    std::vector<typename std::list<LevelCont>::iterator> _LevelStarts;

  public:
    GlobalListCL  (Uint LevelBegin, Uint Size=0)
        : _LevelBegin(LevelBegin), _Data(Size), _LevelStarts(Size)
        { Uint i=0;
          for (iterator It(_Data.begin()); It!=_Data.end(); ++It) _LevelStarts[i++] = It; }
    // standard dtor

    Uint GetLevelBegin   () const { return _LevelBegin; }
    // One greater than the last accessible Level!
    Uint GetLevelEnd     () const { return _LevelBegin+_LevelStarts.size(); }
    Uint GetSize         () const { return _LevelStarts.size(); }
    Uint GetFullSize     () const;
    bool IsEmpty         () const { return _LevelStarts.empty(); }

    const_iterator begin () const { return _Data.begin(); }
          iterator begin ()       { return _Data.begin(); }
    const_iterator end   () const { return _Data.end(); }
          iterator end   ()       { return _Data.end(); }

    void AppendLevel     () { _Data.push_back(std::list<T>()); _LevelStarts.push_back(--_Data.end()); }
    void RemoveLastLevel ()
        { Assert(_Data.back().empty(), DROPSErrCL("LevelList: RemoveLevel: back() not empty"), DebugContainerC);
          _Data.pop_back(); _LevelStarts.pop_back(); }

    const LevelCont&  operator [] (Uint Level) const { return *_LevelStarts[Level-_LevelBegin]; }
    LevelCont&        operator [] (Uint Level)       { return *_LevelStarts[Level-_LevelBegin]; }
};

template <class T>
Uint GlobalListCL<T>::GetFullSize () const
{
    Uint fullsize=0;

    for (const_iterator it(begin()); it!=end(); ++it) fullsize+=it->size();
    return fullsize;
}


} // end of namespace DROPS

#endif
