#include "num/spmat.h"
#include "misc/utils.h"

#include <cmath>

using namespace DROPS;
using namespace std;


template <typename T>
class MyVectorBaseCL: public std::valarray<T>
{
public:
    typedef T value_type;
    typedef std::valarray<T> base_;

    // ctors
    MyVectorBaseCL ()                                 : base_ ()     {}
#ifdef VALARRAY_BUG
    MyVectorBaseCL (size_t s)                         : base_(T(),s) {}
#else
    MyVectorBaseCL (size_t s)                         : base_(s)     {}
#endif
    MyVectorBaseCL (T c, size_t s)                    : base_(c, s)  {}
    MyVectorBaseCL (const T* tp, size_t s)            : base_(tp, s) {}

    template <class V>
    MyVectorBaseCL (V expt)                           : base_(expt)  {}

    // member functions
    T norm2     () const { return (*this)*(*this); }
    T norm      () const { return std::sqrt(norm2()); }
    T supnorm   () const { return std::abs(_va).max(); }
};

typedef MyVectorBaseCL<double> MyVectorCL;

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
  inline double
  TestVector( const VT& v1, const VT& v2, const VT& v3, double a)
{
    VT t1= v1 + 3.0* v2 + a*v3;
    t1/=2.0;
    return (t1+v2)[0];
}


int main()
{   TimerCL time;
    typedef valarray<double> VAT;
    VAT x(VecSize), y(VecSize), z(VecSize);
    double a= -0.2;
    time.Start();
    InitVector( x);InitVector( y);InitVector( z);
    time.Stop();
    cout << "Valarray Init: " << time.GetTime() << endl;
    time.Reset();
    cout << (x+y).max() << endl;

    time.Start();
    double ret= TestVector( x, y, z, a);
    time.Stop();
    cout << "Valarray Calc: " << time.GetTime() << endl;
    time.Reset();


    VectorCL vx(VecSize), vy(VecSize), vz(VecSize);
    time.Start();
    InitVector( vx);InitVector( vy);InitVector( vz);
    time.Stop();
    cout << "VectorCL Init: " << time.GetTime() << endl;
    time.Reset();
    cout << (vx+vy).max() << endl;

    time.Start();
    double vret= TestVector( vx, vy, vz, a);
    time.Stop();
    cout << "VectorCL Calc: " << time.GetTime() << endl;
    time.Reset();

    MyVectorCL myx(VecSize), myy(VecSize), myz(VecSize);
    time.Start();
    InitVector( myx);InitVector( myy);InitVector( myz);
    time.Stop();
    cout << "MyVectorCL Init: " << time.GetTime() << endl;
    time.Reset();

    time.Start();
    double myret= TestVector( myx, myy, myz, a);
    time.Stop();
    cout << "MyVectorCL Calc: " << time.GetTime() << endl;
    time.Reset();
    return static_cast<int>( ret) + static_cast<int>( vret) + static_cast<int>( myret);
}
