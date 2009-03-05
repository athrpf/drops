//**************************************************************************
// File:    quadOut.h                                                      *
// Content: solution output on a quadrilateral grid                        *
//          formaly known as "Volker-patch"                                *
// Author:  Joerg Grande, Sven Gross, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************

#ifndef DROPS_QUADOUT_H
#define DROPS_QUADOUT_H

#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include "misc/utils.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/spmat.h"

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace DROPS
{

template<typename T>
/// \brief Class for storing a three-dimensional array
class Matrix3DCL
{
  private:
    VectorBaseCL<T> val_;       // values
    Uint            dim_[3];    // number of values in each direction

  public:
    Matrix3DCL() :val_(1) {}
    Matrix3DCL(const Matrix3DCL& mat) : val_(mat.val_)
      {
          dim_[0]=mat.dim_[0]; dim_[1]=mat.dim_[1]; dim_[2]=mat.dim_[2];
      }
    Matrix3DCL(Uint nx, Uint ny, Uint nz)
      {
          resize(nx,ny,nz);
      }
    /// \brief resize
    void resize(Uint nx, Uint ny, Uint nz)
      {
          dim_[0]=nx; dim_[1]=ny; dim_[2]=nz;
          val_.resize(nx*ny*nz);
      }
    /// \brief Get number of elements
    Uint size() const
      {
          return val_.size();
      }
    /// \brief Get number of values in each direction
    const Uint* dim() const
      {
          return dim_;
      }
    /// \brief Writing access
    T& operator() (Uint i, Uint j, Uint k)
      {
          Assert(i<dim_[0] && j<dim_[1] && k<dim_[k], DROPSErrCL("Matrix3DCL::operator(): wrong index"), DebugContainerC);
          return val_[k*dim_[0]*dim_[1] + j*dim_[0] +i];
      }
    /// \brief Reading access
    const T& operator() (Uint i, Uint j, Uint k) const
      {
          Assert(i<dim_[0] && j<dim_[1] && k<dim_[k], DROPSErrCL("Matrix3DCL::operator() const: wrong index"), DebugContainerC);
          return val_[k*dim_[0]*dim_[1] + j*dim_[0] +i];
      }
    /// \brief Get constant pointer to all elements
    const T* GetVals() const
      {
          return Addr(val_);
      }
    /// \brief Get pointer to all elements
    T* GetVals()
      {
          return Addr(val_);
      }
};

/// \brief Integrate along one axis to get integrated values on one plane
template<typename T>
void Integrate(Matrix3DCL<T>& output, const T& dx, Uint integration_axis=2, Uint result_plane=0)
/// This function takes the the output and sums all elements along the given
/// dimension up and scales it by dx.
/// \param output The values on the quadrilateral grid of velocities,
///               pressure, ... The result is contained in this matrix
/// \param dx Scaling factor
/// \param integration_axis sum up elements along this axis
/// \param result_plane plane in output, where to store the result
{
    switch (integration_axis)
    {
        // integrate along the x axis
        case 0:
            for (Uint j=0; j<output.dim()[1]; ++j){
                for (Uint k=0; k<output.dim()[2]; ++k){
                    for (Uint i=0; i<output.dim()[0]; ++i)
                        if (i!=result_plane)
                          output(result_plane,j,k) += output(i,j,k);
                    output(result_plane,j,k) *= dx;
                }
            }
            break;
        // integrate along the y axis
        case 1:
            for (Uint i=0; i<output.dim()[0]; ++i){
                for (Uint k=0; k<output.dim()[2]; ++k){
                    for (Uint j=0; j<output.dim()[1]; ++j)
                        if (j!=result_plane)
                          output(i,result_plane,k) += output(i,j,k);
                    output(i,result_plane,k) *= dx;
                }
            }
            break;
        // integrate along the z axis
        case 2:
            for (Uint i=0; i<output.dim()[0]; ++i){
                for (Uint j=0; j<output.dim()[1]; ++j){
                    for (Uint k=0; k<output.dim()[2]; ++k)
                      if (k!=result_plane)
                        output(i,j,result_plane) += output(i,j,k);
                    output(i,j,result_plane) *= dx;
                }
            }
            break;
        default:
            throw DROPSErrCL("Integrate: integration_axis must be 0, 1 or 2");
    };
}

/// \brief Write out norms of the vector-values of an array onto an output-stream
template<typename T>
inline void WriteNormValsOnZPlane(std::ostream& os, T* vals, Uint nx, Uint ny, Uint nz,
                                  Uint format=12, T scale_factor=1.)
/// This function writes out values of an array, that has been created by the
/// Matrix3DCL.
/// \param os output stream
/// \param vals values in an array
/// \param nx number of values in x direction
/// \param ny number of values in y direction
/// \param nz number of values in z direction
/// \param format number of digits to encode the value (default =12)
/// \param scale_factor Scaling factor for the output (default =1.)
{
    Comment("DIMENSIONS IN NORM BERECHNUNG:\n", DebugOutPutC);
    T accu;
    Uint i=0;
    Uint j=0;
    Uint k=0;
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                accu=0;
                for (Uint item=0; item<3; ++item){
                    accu+=vals[3*(k*nx*ny + j*nx + i) + item]*vals[3*(k*nx*ny + j*nx + i) + item];
                }
                os << std::setw(format) << std::sqrt(accu)*scale_factor <<" ";
            }
        }
    }
    Comment(i<<" x "<< j<<" x "<< k<<" \n ", DebugOutPutC);
}

/// \brief Write out the values of a field of scalar values stored in an array
template<typename T>
inline void WriteValsOnZPlane(std::ostream& os, T* vals, Uint nx, Uint ny, __UNUSED__ Uint nz, Uint result_plane=0,
                              Uint format=12, Uint numperitem=1, int single_item=-1, T scale_factor=1)
/// This function writes out values of an array, that has been created by the
/// Matrix3DCL.
/// \param os output stream
/// \param vals values in an array
/// \param nx number of values in x direction
/// \param ny number of values in y direction
/// \param nz number of values in z direction
/// \param result_plane
/// \param format number of digits to encode the value (default =12)
/// \param numperitem number of numbers per value (3 for velocity) (default= 1)
/// \param single_item if specified only the single_items's number of the value
///                    is printed out (default=-1 all values are printed)
/// \param scale_factor Scaling factor for the output (default =1.)
{
    if (single_item>=(int)numperitem)
        throw DROPSErrCL("WriteVals: Wrong numperitem and single_item has been specified");
    Uint k=result_plane;
    for (Uint j=0; j<ny; ++j) {
        for (Uint i=0; i<nx; ++i){
            if (single_item<0){
                for (Uint item=0; item<numperitem; ++item){
                    os << std::setw(format) << vals[numperitem*(k*nx*ny + j*nx + i) + item]*scale_factor << ' ';
                }
            }
            else{
                os << std::setw(format) <<vals[numperitem*(k*nx*ny + j*nx + i) +single_item]*scale_factor<<' ';
            }
        }
        os << '\n';
    }
    os << '\n';
}

/// \brief Write out the values of a field of scalar values stored in an array
template<typename T>
inline void WriteVals(std::ostream& os, T* vals, Uint nx, Uint ny, Uint nz, Uint format=12,
                      Uint numperitem=1, int single_item=-1, T scale_factor=1)
/// This function writes out values of an array, that has been created by the
/// Matrix3DCL.
/// \param os output stream
/// \param vals values in an array
/// \param nx number of values in x direction
/// \param ny number of values in y direction
/// \param nz number of values in z direction
/// \param format number of digits to encode the value (default =12)
/// \param numperitem number of numbers per value (3 for velocity) (default= 1)
/// \param single_item if specified only one number per value is printed out
///                   (default=-1 all values are printed)
/// \param scale_factor Scaling factor for the output (default =1.)
{
    if (single_item>=(int)numperitem)
        throw DROPSErrCL("WriteVals: Wrong numperitem and single_item has been specified");
    for (Uint k=0; k<nz; ++k) {
        os << "\n# Layer " << k << '\n';
        for (Uint j=0; j<ny; ++j) {
            for (Uint i=0; i<nx; ++i){
                if (single_item<0){
                    for (Uint item=0; item<numperitem; ++item)
                        os << std::setw(format) << vals[numperitem*(k*nx*ny + j*nx + i) + item]*scale_factor << ' ';
                }
                else{
                    os << std::setw(format) <<vals[numperitem*(k*nx*ny + j*nx + i) +single_item]*scale_factor<<' ';
                }
            }
            os << '\n';
        }
    }
    os << '\n';
}

/// \brief Transform a Matrix3DCL of 3D-points to an array
inline double* Mat3DPoints2Double(const Matrix3DCL<Point3DCL>& mat)
/** This function takes a Matrix3DCL of 3D-points and returns a pointer to an
    array of all values*/
{
    double * result= new double[3*mat.size()];
    const Point3DCL *p=mat.GetVals();
    for (Uint i=0; i<mat.size(); ++i, ++p)
         for (int j=0; j<3; ++j)
       result[3*i+j] = (*p)[j];
    return result;
}

/// \brief Class for writing a solution on a quadrilateral grid
class QuadOutCL
{
  private:
    const MultiGridCL& MG_;
    const IdxDescCL*   idx_;
    int                gridpts_[3];  // grid points in x-, y-, z-direction
    Point3DCL          h_;           // step size in each direction
    Point3DCL          bary_;        // barycenter of the grid
    Point3DCL          offset_;
    Point3DCL          rot_;         // rotation in each direction
    SMatrixCL<3,3>     rotmat_;      // rotation matrix determined by rot_
#ifdef _PAR
    Matrix3DCL<int>   *geom_info_;   // this matrix contains information, how many processors write out a single point of the quadrilateral gird
    size_t             geom_version_;// version of the multigrid (for geom_info_)
    Matrix3DCL<LocationCL> locations_;// store information, where the quadrilateral can be found
    ProcCL::OperationT reduceOp_;     // operation for user defined reduce
#endif

  private:
    void CheckFile( const std::ofstream&) const;

#ifdef _PAR
    template <typename T>
    static inline void PerformReduce(const void* invec, void* inoutvec, int len);
    static inline void ReduceCXX(const void* invec, void* inoutvec, int len, const ProcCL::DatatypeT& datatype);
    static inline void ReduceC(void* invec, void* inoutvec, int* len, ProcCL::DatatypeT* datatype);
#endif

  public:
    QuadOutCL( const MultiGridCL& mg, const IdxDescCL* idx)
      : MG_( mg), idx_( idx), geom_info_(0), geom_version_(size_t(-1))
    /** \note: idx must refer to a numbering of the verts and edges of a certain
        triangulation, i.e. use LevelsetP2CL::CreateNumbering before
        constructing an QuadOutCL-object. */
      { Init(1,1,1,Point3DCL(),Point3DCL(),Point3DCL()); }
    QuadOutCL( const MultiGridCL& mg, const IdxDescCL* idx,
               int nx, int ny, int nz, const Point3DCL& h,
               const Point3DCL& bary, const Point3DCL& rot)
      : MG_( mg), idx_( idx), geom_info_(0), geom_version_(size_t(-1))
      { Init(nx, ny, nz, h, bary, rot);}
    ~QuadOutCL() {
        if (geom_info_) delete geom_info_;
#ifdef _PAR
        ProcCL::FreeOp(reduceOp_);
#endif
    }

    /// \brief Store geometric information about the quadrilateral grid
    void Init(int nx, int ny, int nz, const Point3DCL& h,
              const Point3DCL& bary, const Point3DCL& rot);

    /// \brief Write information about quadrilateral grid in a file
    void putLog    ( std::string);
    /// \brief Write geometry of the quadrilateral grid in a file
    void putGeom   ( std::string);
    /// \brief Write a scalar valued finite element function on a quadrilateral grid in a file
    template<class DiscScalT>
    void putScalar ( std::string, const DiscScalT&);
    /// \brief Write a vector valued finite element function on a quadrilateral grid in a file
    template<typename DiscVecT, typename DescScalT>
    void putVector (  std::string, std::string, std::string, const DiscVecT&, const DescScalT* lset);
    /// \brief Write a vector valued finite element function on a quadrilateral grid in a file
    template<typename DiscVecT>
    void putVector (  std::string, std::string, std::string, const DiscVecT&);

};


inline void QuadOutCL::Init(int nx, int ny, int nz, const Point3DCL& h,
                            const Point3DCL& bary, const Point3DCL& rot)
{
    gridpts_[0]= nx;
    gridpts_[1]= ny,
    gridpts_[2]= nz;
    offset_[0]= (nx-1)/2.;
    offset_[1]= (ny-1)/2.;
    offset_[2]= (nz-1)/2.;
    h_= h;
    offset_*= h_;
    bary_= bary;
    rot_= rot;

    SMatrixCL<3,3> RX, RY, RZ;

    RX(0,0)= RY(1,1)= RZ(2,2)= 1.;
    RX(1,1)=   RX(2,2)= std::cos(rot_[0]);
    RX(2,1)= -(RX(1,2)= std::sin(rot_[0]));
    RY(0,0)=   RY(2,2)= std::cos(rot_[1]);
    RY(2,0)= -(RY(0,2)= std::sin(rot_[1]));
    RZ(0,0)=   RZ(1,1)= std::cos(rot_[2]);
    RZ(1,0)= -(RZ(0,1)= std::sin(rot_[2]));

    rotmat_= RZ * RY * RX;

    locations_.resize(gridpts_[0], gridpts_[1], gridpts_[2]);
#ifdef _PAR
# ifdef _MPICXX_INTERFACE
    ProcCL::InitOp(reduceOp_, &ReduceCXX, true);
# else
    ProcCL::InitOp(reduceOp_, &ReduceC, true);
#endif
#endif
}


inline void QuadOutCL::CheckFile( const std::ofstream& os) const
{
    if (!os) throw DROPSErrCL( "QuadOutCL: error while opening file!");
}


#ifdef _PAR
template <typename T>
inline void QuadOutCL::PerformReduce(const void* invec, void* inoutvec, int len)
{
    T* in   = (T*) invec;
    T* inout= (T*) inoutvec;
    for (int i=0; i<len; ++i){
    	const T abs_inout= inout[i]<0 ? -inout[i] : inout[i];
    	const T abs_in   = in[i]<0    ? -in[i]    : in[i];
        inout[i]= abs_inout>abs_in ? inout[i] : in[i];
	}
}

inline void QuadOutCL::ReduceC(void* invec, void* inoutvec, int* len, ProcCL::DatatypeT* datatype)
{
    if ( *datatype==ProcCL::MPI_TT<double>::dtype ){
        PerformReduce<double>(invec, inoutvec, *len);
    }
    else if ( *datatype==ProcCL::MPI_TT<int>::dtype  ){
        PerformReduce<int>(invec, inoutvec, *len);
    }
}

inline void QuadOutCL::ReduceCXX(const void* invec, void* inoutvec, int len, const ProcCL::DatatypeT& datatype)
{
    if ( datatype==ProcCL::MPI_TT<double>::dtype ){
        PerformReduce<double>(invec, inoutvec, len);
    }
    else if ( datatype==ProcCL::MPI_TT<int>::dtype  ){
        PerformReduce<int>(invec, inoutvec, len);
    }
}
#endif

inline void QuadOutCL::putLog( std::string fileName)
{
    std::ofstream os( fileName.c_str());
    CheckFile( os);

    os << "Gridpoints: "
       << gridpts_[0] << ' ' << gridpts_[1] << ' ' << gridpts_[2]
       << "\nStepsize: "
       << h_[0] << ' ' << h_[1] << ' ' << h_[2]
       << "\nBarycenter: "
       << bary_[0] << ' ' << bary_[1] << ' ' << bary_[2]
       << "\nRotation angles: "
       << rot_[0] << ' ' << rot_[1] << ' ' << rot_[2]
       << '\n';
}


inline void QuadOutCL::putGeom( std::string fileName)
/** This function writes zeros and ones in the file given by the parameter
    fileName if the corresponding point in the quadrilateral grid is in- or
    outside the computational domain.<p>
    In difference to the sequential version of DROPS, these information are
    gathered by each processor but collected and written out by the master
    processor.<p>
    In the parallel version, this function has to be called before writing out
    vectorial or scalar valued functions. It may happen, that a single point
    of the quadrilateral is written out by two processors. So the value of the
    function is divided by the number of processors, that writes out the value.
\note If DROPS is compiled with OpenMP support, this function is parallelized*/
{
#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif
    double duration=0;

    const Uint lvl= idx_->TriangLevel();
#ifdef _PAR
    Matrix3DCL<int> output(gridpts_[0], gridpts_[1], gridpts_[2]);
#endif

    std::ofstream *os=0;
    IF_MASTER
    {
        os = new std::ofstream( fileName.c_str());
        CheckFile( *os);
        os->flags(std::ios_base::scientific);
        os->precision(7);
        (*os) << "# DROPS data file, scalar variable\n# grid dimensions:\n"
              << gridpts_[0] << " " << gridpts_[1] << " " << gridpts_[2] << "\n";
    }

    // In sequential mode of DROPS: write out zeroes and ones into a file,
    // in parallel mode of DROPS: Collect information in a Matrix3DCL<int>
#pragma omp parallel
{
#ifdef _OPENMP
    int num_threads=omp_get_num_threads();
#else
    int num_threads=1;
#endif
    IF_MASTER{
#pragma omp master
        std::cout << "   * Using "<<num_threads<<" thread(s) to get geometry ..." << std::endl;
    }
#pragma omp for schedule(dynamic)
    for (int k=0; k<gridpts_[2]; ++k)
    {
#ifndef _PAR
        (*os) << "\n# Layer " << k << '\n';
#endif
        for (int j=0; j<gridpts_[1]; ++j)
        {
            for (int i=0; i<gridpts_[0]; ++i)
            {
                Point3DCL p= MakePoint3D( i*h_[0], j*h_[1], k*h_[2]);
                p= rotmat_ * (p - offset_) + bary_;
#ifndef _PAR
                LocationCL loc;
                LocatorCL::Locate( loc, MG_, lvl, p);
                (*os) << std::setw(2) << loc.IsValid(lvl);
#else
                LocatorCL::Locate( locations_(i,j,k), MG_, lvl, p);
                output(i,j,k)= locations_(i,j,k).IsValid(lvl);
#endif
            }
#ifndef _PAR
            (*os) << '\n';
#endif
        }
    }
}       // end of openmp

#ifdef _PAR
    // Collect information on master processor, which writes values in the given file
    geom_info_   = new Matrix3DCL<int>(output);
    geom_version_= MG_.GetVersion();
    ProcCL::GlobalOp(output.GetVals(), geom_info_->GetVals(), output.size(), -1, reduceOp_);
    if (ProcCL::IamMaster())
        WriteVals(*os, geom_info_->GetVals(), gridpts_[0], gridpts_[1], gridpts_[2], 2);
#endif
    if (os) delete os;
    timer.Stop(); duration=timer.GetTime();
    IF_MASTER
        std::cerr << "   * took "<<duration<<" sec."<<std::endl;
}


template<class DiscScalT>
void QuadOutCL::putScalar( std::string fileName, const DiscScalT& v)
/** This function writes out a scalar valued finite element function v on a
    quadrilateral grid.<p>
    In difference to the sequential version of DROPS, these information are
    gathered by each processor but collected and written out by the master
    processor.
\note If DROPS is compiled with OpenMP support, this function is parallelized
\pre The function putGeom has to be called for the same multigrid.*/
{
    if (!geom_info_ || geom_version_!=MG_.GetVersion())
        throw DROPSErrCL("QuadOutCL::putScalar: Call putGeom before calling putScalar");
    // This function works like puGeom. For more comments look there
    const Uint lvl= idx_->TriangLevel;
#ifdef _PAR
    Matrix3DCL<double> output(gridpts_[0], gridpts_[1], gridpts_[2]);
#endif
    std::ofstream *os=0;
    IF_MASTER
    {
        os = new std::ofstream( fileName.c_str());
        CheckFile( *os);
        os->flags(std::ios_base::scientific);
        os->precision(7);
        (*os) << "# DROPS data file, scalar variable\n# grid dimensions:\n"
              << gridpts_[0] << " " << gridpts_[1] << " " << gridpts_[2] << "\n";
    }

    for (int k=0; k<gridpts_[2]; ++k)
    {
#ifndef _PAR
        (*os) << "\n# Layer " << k << '\n';
#endif
        for (int j=0; j<gridpts_[1]; ++j)
        {
            for (int i=0; i<gridpts_[0]; ++i)
            {
#ifndef _PAR
                Point3DCL p= MakePoint3D( i*h_[0], j*h_[1], k*h_[2]);
                p= rotmat_ * (p - offset_) + bary_;

                LocationCL loc;
                LocatorCL::Locate( loc, MG_, lvl, p);
#else
                LocationCL &loc=locations_(i,j,k);
#endif
                const double val = loc.IsValid(lvl)
                                 ? v.val( loc.GetTetra(), loc.GetBaryCoord())
                                 : 0;
#ifndef _PAR
                (*os) << std::setw(12) << val << ' ';
#else
                output(i,j,k) = val / std::max(1,(*geom_info_)(i,j,k));
#endif
            }
#ifndef _PAR
            (*os) << '\n';
#endif
        }
    }

#ifdef _PAR
    Matrix3DCL<double> *globalsol=0;
    if (ProcCL::IamMaster())
        globalsol = new Matrix3DCL<double>(output);
    ProcCL::GlobalOp(output.GetVals(),
                     ProcCL::IamMaster()? globalsol->GetVals() : 0,
                     output.size(), ProcCL::Master(), reduceOp_);
    if (ProcCL::IamMaster())
        WriteVals(*os, globalsol->GetVals(), gridpts_[0], gridpts_[1], gridpts_[2], 12);
    if (globalsol) delete globalsol;
#endif
    if (os) delete os;
}


template<typename DiscVecT, typename DescScalT>
void QuadOutCL::putVector(std::string fileName, std::string fileNameY, std::string fileNameZ, const DiscVecT& v, const DescScalT* lset)
/** Note: If DROPS is compiled with OpenMP support, this function is
    parallelized. If the parameter lset points to a level-set function, only
    values are written out, where the level-set values are smaller or equal zero
\pre The function putGeom has to be called for the same multigrid.*/
{
    if (!geom_info_ || geom_version_!=MG_.GetVersion())
        throw DROPSErrCL("QuadOutCL::putVector: Call putGeom before calling putVector");

#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif
    double duration=0;

    Comment("PUT VECTOR\n", DebugOutPutC);
    const Uint lvl= idx_->TriangLevel();
    Matrix3DCL<Point3DCL> output(gridpts_[0], gridpts_[1], gridpts_[2]);

    std::ofstream *os=0;
    std::ofstream *osY=0;
    std::ofstream *osZ=0;

    IF_MASTER
    {
        os = new std::ofstream( fileName.c_str());
        CheckFile( *os);
        os->flags(std::ios_base::scientific);
        os->precision(7);

        osY = new std::ofstream( fileNameY.c_str());
        CheckFile( *osY);
        osY->flags(std::ios_base::scientific);
        osY->precision(7);

        osZ = new std::ofstream( fileNameZ.c_str());
        CheckFile( *osZ);
        osZ->flags(std::ios_base::scientific);
        osZ->precision(7);
    }

    Comment("Fetch the points\n", DebugOutPutC);

#pragma omp parallel
{
#ifdef _OPENMP
    int num_threads=omp_get_num_threads();
#else
    int num_threads=1;
#endif

    IF_MASTER{
#pragma omp master
        std::cout << "   * Using "<<num_threads<<" thread(s) to fetch points ..." << std::endl;
    }
#pragma omp for
    for (int k=0; k<gridpts_[2]; ++k)
    {
        for (int j=0; j<gridpts_[1]; ++j)
        {
            for (int i=0; i<gridpts_[0]; ++i)
            {
#ifndef _PAR
                Point3DCL p= MakePoint3D( i*h_[0], j*h_[1], k*h_[2]);
                p= rotmat_ * (p - offset_) + bary_;

                LocationCL loc;
                LocatorCL::Locate( loc, MG_, lvl, p);
#else
                LocationCL &loc=locations_(i,j,k);
#endif
                Point3DCL val;
                if (lset!=0)
                    val = loc.IsValid(lvl) && lset->val(loc.GetTetra(), loc.GetBaryCoord())<=0
                        ? v.val( loc.GetTetra(), loc.GetBaryCoord())
                        : Point3DCL();
                else
                    val = loc.IsValid(lvl)
                            ? v.val( loc.GetTetra(), loc.GetBaryCoord())
                    : Point3DCL();

                output(i,j,k) = val / std::max(1,(*geom_info_)(i,j,k));
            }
        }
    }
}       // end of openmp
    Comment("putVector: Start_integrate  ...  ", DebugOutPutC);
    Integrate(output, h_);

    Comment("DONE\n putVector: output to array of double  ...  ", DebugOutPutC);
    double *localsol  = Mat3DPoints2Double(output);

    double *globalsol=0;
#ifndef _PAR
    globalsol = localsol;
#else
    if (ProcCL::IamMaster())
        globalsol = new double[3*output.size()];
    ProcCL::GlobalOp(localsol,
                     ProcCL::IamMaster()? globalsol : 0,
                     3*output.size(), ProcCL::Master(), reduceOp_);
#endif

    IF_MASTER
    {
        Comment("DONE\n putVector: writenorm   ...  ", DebugOutPutC);
        WriteNormValsOnZPlane(*os, globalsol, gridpts_[0], gridpts_[1], gridpts_[2], 12, 1.);

        Comment("DONE\n putVector: write Z plane with y value of the velocity ...  ", DebugOutPutC);
        WriteValsOnZPlane(*osY, globalsol, gridpts_[0], gridpts_[1], gridpts_[2],  0, 12, 3, 1, 1.);

        Comment("DONE\n putVector: write Z plane with z value of the velocity ...  ", DebugOutPutC);
        WriteValsOnZPlane(*osZ, globalsol, gridpts_[0], gridpts_[1], gridpts_[2],  0, 12, 3, 2, 1.);
        Comment("DONE\n", DebugOutPutC);
    }

    // Clean up
    if (globalsol) delete globalsol;
    if (os) delete os;
    if (osY) delete osY;
    if (osZ) delete osZ;
    timer.Stop(); duration=timer.GetTime();
    IF_MASTER
        std::cerr << "   * took "<<duration<<" sec."<<std::endl;

}

template<typename DiscVecT>
void QuadOutCL::putVector(std::string fileName, std::string fileNameY, std::string fileNameZ, const DiscVecT& v)
/** Note: If DROPS is compiled with OpenMP support, this function is
    parallelized. If the parameter lset points to a level-set function, only
    values are written out, where the level-set values are smaller or equal zero
\pre The function putGeom has to be called for the same multigrid.*/
{
    if (!geom_info_ || geom_version_!=MG_.GetVersion())
        throw DROPSErrCL("QuadOutCL::putVector: Call putGeom before calling putVector");
#ifndef _PAR
    TimerCL timer;
#else
    ParTimerCL timer;
#endif
    double duration=0;

    Comment("PUT VECTOR\n", DebugOutPutC);
    const Uint lvl= idx_->TriangLevel();
    Matrix3DCL<Point3DCL> output(gridpts_[0], gridpts_[1], gridpts_[2]);

    std::ofstream *os=0;
    std::ofstream *osY=0;
    std::ofstream *osZ=0;

    IF_MASTER
    {
        os = new std::ofstream( fileName.c_str());
        CheckFile( *os);
        os->flags(std::ios_base::scientific);
        os->precision(7);

        osY = new std::ofstream( fileNameY.c_str());
        CheckFile( *osY);
        osY->flags(std::ios_base::scientific);
        osY->precision(7);

        osZ = new std::ofstream( fileNameZ.c_str());
        CheckFile( *osZ);
        osZ->flags(std::ios_base::scientific);
        osZ->precision(7);
    }

    Comment("Fetch the points\n", DebugOutPutC);

#pragma omp parallel
{
#ifdef _OPENMP
    int num_threads=omp_get_num_threads();
#else
    int num_threads=1;
#endif

    IF_MASTER{
#pragma omp master
        std::cerr << "   * Using "<<num_threads<<" thread(s) to fetch points ..." << std::endl;
    }
#pragma omp for
    for (int k=0; k<gridpts_[2]; ++k)
    {
        for (int j=0; j<gridpts_[1]; ++j)
        {
            for (int i=0; i<gridpts_[0]; ++i)
            {
#ifndef _PAR
                Point3DCL p= MakePoint3D( i*h_[0], j*h_[1], k*h_[2]);
                p= rotmat_ * (p - offset_) + bary_;

                LocationCL loc;
                LocatorCL::Locate( loc, MG_, lvl, p);
#else
                LocationCL &loc=locations_(i,j,k);
#endif
                Point3DCL val;
                val = loc.IsValid(lvl)
                        ? v.val( loc.GetTetra(), loc.GetBaryCoord())
                        : Point3DCL();

                output(i,j,k) = val / std::max(1,(*geom_info_)(i,j,k));
            }
        }
    }
}       // end of openmp

    Comment("putVector: Start_integrate  ...  ", DebugOutPutC);
    Integrate(output, h_);

    Comment("DONE\n putVector: output to array of double  ...  ", DebugOutPutC);
    double *localsol  = Mat3DPoints2Double(output);

    double *globalsol=0;
#ifndef _PAR
    globalsol = localsol;
#else
    if (ProcCL::IamMaster())
        globalsol = new double[3*output.size()];
    ProcCL::GlobalOp(localsol,
                     ProcCL::IamMaster()? globalsol : 0,
                     3*output.size(), ProcCL::Master(), reduceOp_);
#endif

    IF_MASTER
    {
        Comment("DONE\n putVector: writenorm   ...  ", DebugOutPutC);
        WriteNormValsOnZPlane(*os, globalsol, gridpts_[0], gridpts_[1], gridpts_[2], 12, 1.);

        Comment("DONE\n putVector: write Z plane with y value of the velocity ...  ", DebugOutPutC);
        WriteValsOnZPlane(*osY, globalsol, gridpts_[0], gridpts_[1], gridpts_[2],  0, 12, 3, 1, 1.);

        Comment("DONE\n putVector: write Z plane with z value of the velocity ...  ", DebugOutPutC);
        WriteValsOnZPlane(*osZ, globalsol, gridpts_[0], gridpts_[1], gridpts_[2],  0, 12, 3, 2, 1.);
        Comment("DONE\n", DebugOutPutC);
    }

    // Clean up
    if (globalsol) delete globalsol;
    if (os) delete os;
    if (osY) delete osY;
    if (osZ) delete osZ;
    timer.Stop(); duration=timer.GetTime();
    IF_MASTER
        std::cerr << "   * took "<<duration<<" sec."<<std::endl;
}

} // end of namespace DROPS

#endif
