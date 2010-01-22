/// \file vectest.cpp
/// \brief compares old and new VectorCL-implementation
/// \author LNM RWTH Aachen: Joerg Grande, Volker Reichelt; SC RWTH Aachen:

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

#include "num/spmat.h"
#include "misc/utils.h"

#include <cmath>

#ifdef __SUNPRO_CC
#  include <stdlib.h>
#endif

using namespace DROPS;
using namespace std;


template <typename T>
class OldVectorBaseCL
{
private:
    std::valarray<T> _va;

    void OldAssertDim (__UNUSED__ const std::valarray<T>& va, __UNUSED__ const char msg[]) const
      { Assert(_va.size()==va.size(), msg, DebugNumericC); }

public:
    typedef T value_type;

    // ctors
    OldVectorBaseCL ()                                 : _va()      {}
    OldVectorBaseCL (size_t s)                         : _va(s)     {}
    OldVectorBaseCL (T c, size_t s)                    : _va(c, s)  {}
    OldVectorBaseCL (const T* tp, size_t s)            : _va(tp, s) {}
    OldVectorBaseCL (const std::valarray<T>& va)       : _va(va)    {}
    OldVectorBaseCL (const OldVectorBaseCL& v)            : _va(v._va) {}

    OldVectorBaseCL (const std::slice_array<T>& sla)   : _va(sla)   {}
    OldVectorBaseCL (const std::gslice_array<T>& gsla) : _va(gsla)  {}
    OldVectorBaseCL (const std::mask_array<T>& ma)     : _va(ma)    {}
    OldVectorBaseCL (const std::indirect_array<T>& ia) : _va(ia)    {}

    void resize (size_t s, T c = T()) { _va.resize(s, c); }

    const std::valarray<T>& raw() const {
        return _va;}
    std::valarray<T>& raw() {
        return _va;}

    // element access
    T  operator [] (size_t s) const
      { Assert(s<size(), "OldVectorBaseCL []: index out of bounds", DebugNumericC); return _va[s]; }
    T& operator [] (size_t s)
      { Assert(s<size(), "OldVectorBaseCL []: index out of bounds", DebugNumericC); return _va[s]; }

    // assignment
    OldVectorBaseCL& operator= (const std::valarray<T>& va)
      { OldAssertDim(va,   "OldVectorBaseCL =: incompatible dimensions"); _va = va;    return *this; }
    OldVectorBaseCL& operator= (const OldVectorBaseCL& v)
      { OldAssertDim(v._va,"OldVectorBaseCL =: incompatible dimensions"); _va = v._va; return *this; }

    OldVectorBaseCL& operator= (T c)                              { _va = c;    return *this; }
    OldVectorBaseCL& operator= (const std::slice_array<T>& sla)   { _va = sla;  return *this; }
    OldVectorBaseCL& operator= (const std::gslice_array<T>& gsla) { _va = gsla; return *this; }
    OldVectorBaseCL& operator= (const std::mask_array<T>& ma)     { _va = ma;   return *this; }
    OldVectorBaseCL& operator= (const std::indirect_array<T>& ia) { _va = ia;   return *this; }

    // computed assignment
    OldVectorBaseCL& operator+= (T c) { _va += c; return *this; }
    OldVectorBaseCL& operator-= (T c) { _va -= c; return *this; }
    OldVectorBaseCL& operator*= (T c) { _va *= c; return *this; }
    OldVectorBaseCL& operator/= (T c) { _va /= c; return *this; }
    OldVectorBaseCL& operator+= (const OldVectorBaseCL& v)
      { OldAssertDim(v._va,"OldVectorBaseCL +=: incompatible dimensions"); _va += v._va; return *this; }
    OldVectorBaseCL& operator-= (const OldVectorBaseCL& v)
      { OldAssertDim(v._va,"OldVectorBaseCL -=: incompatible dimensions"); _va -= v._va; return *this; }
    OldVectorBaseCL& operator*= (const OldVectorBaseCL& v)
      { OldAssertDim(v._va,"OldVectorBaseCL *=: incompatible dimensions"); _va *= v._va; return *this; }
    OldVectorBaseCL& operator/= (const OldVectorBaseCL& v)
      { OldAssertDim(v._va,"OldVectorBaseCL /=: incompatible dimensions"); _va /= v._va; return *this; }

    // unary minus
    OldVectorBaseCL operator- () const { return OldVectorBaseCL(-_va); }

    // member functions
    size_t size () const { return _va.size(); }
    T sum       () const { return _va.sum(); }
    T min       () const { return _va.min(); }
    T max       () const { return _va.max(); }
    T norm2     () const { return (*this)*(*this); }
    T norm      () const { return std::sqrt(norm2()); }
    T supnorm   () const { return std::abs(_va).max(); }

    friend OldVectorBaseCL operator+ (const OldVectorBaseCL& v, const OldVectorBaseCL& w)
      { v.OldAssertDim(w._va,"OldVectorBaseCL + OldVectorBaseCL: incompatible dimensions"); return OldVectorBaseCL(v._va+w._va); }
    friend OldVectorBaseCL operator- (const OldVectorBaseCL& v, const OldVectorBaseCL& w)
      { v.OldAssertDim(w._va,"OldVectorBaseCL - OldVectorBaseCL: incompatible dimensions"); return OldVectorBaseCL(v._va-w._va); }
    friend T              operator* (const OldVectorBaseCL& v, const OldVectorBaseCL& w)
    {   v.OldAssertDim(w._va,"OldVectorBaseCL * OldVectorBaseCL: incompatible dimensions");
        return std::inner_product(Addr(v._va),Addr(v._va)+v.size(),Addr(w._va),T());
    }

    friend OldVectorBaseCL operator* (T c, const OldVectorBaseCL& v) { return OldVectorBaseCL(v._va*c); }
    friend OldVectorBaseCL operator* (const OldVectorBaseCL& v, T c) { return OldVectorBaseCL(v._va*c); }
    friend OldVectorBaseCL operator/ (const OldVectorBaseCL& v, T c) { return OldVectorBaseCL(v._va/c); }

    friend void axpy (T a, const OldVectorBaseCL& x, OldVectorBaseCL& y)
    {
        Assert(x.size()==y.size(), "axpy: incompatible dimensions", DebugNumericC);
        y._va+= a*x._va; // y+= a*x;
    }
    friend void z_xpay (OldVectorBaseCL& z, const OldVectorBaseCL& x, T a, const OldVectorBaseCL& y)
    {
        Assert(z.size()==x.size() && z.size()==y.size(), "z_xpay: incompatible dimensions", DebugNumericC);
        z._va= x._va+a*y._va; // z= x+a*y;
    }
    friend void z_xpaypby2 (OldVectorBaseCL& z, const OldVectorBaseCL& x, T a, const OldVectorBaseCL& y, T b, const OldVectorBaseCL& y2)
    {
        Assert(z.size()==x.size() && z.size()==y.size() && z.size()==y2.size(), "z_xpaypby2: incompatible dimensions", DebugNumericC);
        z._va= x._va+a*y._va+b*y2._va; // z= x+a*y+b*y2;
    }
};

typedef OldVectorBaseCL<double> MyVectorCL;


const int VecSize= 10000000;

template <class VT>
  void
  InitVector( VT& v)
{
    for (Uint i= 0; i<v.size(); ++i) {
        v[i]=  drand48();
    }
}

template <class VT>
  inline VT
  TestVector( const VT& v1, const VT& v2, const VT& v3, double a)
{
    VT t1( v1 + 3.0* v2 + a*v3);
    t1/=2.0;
    return VT( t1+v2);
}


int main()
{   TimerCL time;
    typedef valarray<double> VAT;
    VAT x(VecSize), y(VecSize), z(VecSize);
    cout << x.max() << endl;

    double a= -0.2;
    time.Start();
    InitVector( x);InitVector( y);InitVector( z);
    time.Stop();
    cout << "Valarray Init: " << time.GetTime() << endl;
    time.Reset();

    time.Start();
    VAT ret= TestVector( x, y, z, a);
    time.Stop();
    cout << "Valarray Calc: " << time.GetTime() << endl
         << "         Norm: " << norm( ret) << endl;
    time.Reset();

//    VectorCL vx(VecSize), vy(VecSize), vz(VecSize);
    time.Start();
//    InitVector( vx);InitVector( vy);InitVector( vz);
    VectorCL vx( x), vy( y), vz( z);
    time.Stop();
    cout << "VectorCL Init: " << time.GetTime() << endl;
    time.Reset();

    time.Start();
    VectorCL vret= TestVector( vx, vy, vz, a);
    time.Stop();
    cout << "VectorCL Calc: " << time.GetTime() << endl
         << "         Differenz: " << norm( static_cast<VAT>( vret) - ret) << endl;
    time.Reset();

//    MyVectorCL myx(VecSize), myy(VecSize), myz(VecSize);
    time.Start();
//    InitVector( myx);InitVector( myy);InitVector( myz);
    MyVectorCL myx( x), myy( y), myz( z);
    time.Stop();
    cout << "MyVectorCL Init: " << time.GetTime() << endl;
    time.Reset();

    time.Start();
    MyVectorCL myret= TestVector( myx, myy, myz, a).raw();
    time.Stop();
    cout << "MyVectorCL Calc: " << time.GetTime() << endl
         << "           Differenz: " << norm( myret.raw() - ret) << endl;
    time.Reset();
    return static_cast<int>( ret[0]) + static_cast<int>( vret[0])
        + static_cast<int>( myret[0]);
}
