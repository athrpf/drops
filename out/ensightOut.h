//**************************************************************************
// File:    ensightOut.h                                                   *
// Content: solution output in Ensight6 Case format                        *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************

/// \file
/// \todo (merge) EnsightOutCL differs a lot!

#ifndef DROPS_ENSIGHTOUT_H
#define DROPS_ENSIGHTOUT_H

#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <map>
#include "geom/multigrid.h"
#include "misc/problem.h"

#ifdef _PAR
#  include <deque>
#  include <queue>
#  include "parallel/parallel.h"
#  include "parallel/exchange.h"
#endif

namespace DROPS
{

/// \brief Helper union for writing integers in binary files
union showInt
{
    int i;
    char s[sizeof(int)];
};

/// \brief Helper union for writing floats in binary files
union showFloat
{
    float f;
    char s[sizeof(float)];
};


#ifndef _PAR
class Ensight6VariableCL; //forward declaration

/// \brief Class for writing out results of a simulation in Ensight6 Case format.
///
/// Register subclasses of Ensight6VariableCL to output the geometry and scalar/vector-valued functions
/// in Ensight6-Case format.
class Ensight6OutCL
{
  private:
    char               decDigits_; ///< Number of digits in the decimal representation of numsteps_
    int                timestep_;  ///< Current timestep
    int                numsteps_;  ///< Total number of timesteps (required for the names of transient output files)
    double             time_;      ///< Time of current timestep
    std::ostringstream geomdesc_,  ///< Geometry-section of the Case-file
                       vardesc_,   ///< Variable-section of the Case-file
                       timestr_;   ///< Time-section of the Case-file
    std::string        casefile_;  ///< Name of the Case-file
    const bool         binary_;    ///< type of output
    bool               timedep_;   ///< true, if there are time-dependent variables
    std::map<std::string, Ensight6VariableCL*> vars_;        ///< The variables and geometry stored by varName.

    /// \brief Internal helper
    ///@{
    void OpenFile       (std::ofstream& of, std::string varName); ///< Append timecode for transient output and check the stream
    bool putTime        (double t);                               ///< Advance timestep_, time_, timestr_, if t_> time_ and set time_= t_; returns true, if t_>time_.
    void CheckFile      (const std::ofstream&) const;
    void CommitCaseFile ();                                       ///< (Re)write case file
    ///@}

  public:
    Ensight6OutCL  (std::string casefileName, Uint numsteps= 0, bool binary= true);
    ~Ensight6OutCL ();

    /// \brief Register a variable or the geometry for output with Write().
    ///
    /// The class takes ownership of the objects, i. e. it destroys them with delete in its destructor.
    void Register (Ensight6VariableCL& var);
    /// \brief Write the registered Ensight6-variables
    ///
    /// For t==0, write all registered objects to their files;if t>0 and t has increased with regard
    /// to the last call, write all time-dependent objects. Only the first of multiple calls with identical t has an effect.
    void Write (double t= 0.);

    /// \brief Append the current timestep-value with the required number of leading '0' to str.
    void AppendTimecode(std::string& str) const;

    /// \brief Interface for classes that implement the Ensight6VariableCL-interface, i.e. output of specific varibles; should probably be private.
    ///@{
    /// \brief Describe a geometry model
    void DescribeGeom (std::string geoName);
    /// \brief Describe a finite element function
    void DescribeVariable (std::string varName, bool isscalar);

    /// \brief Write the geometry into a file
    void putGeom   (MultiGridCL& mg, int lvl, std::string geoName);
    /// \brief Write a scalar value finite element function into a file
    template<class DiscScalT>
    void putScalar (const DiscScalT& v, std::string varName);
    /// \brief Write a vector value finite element function into a file
    template<class DiscVecT>
    void putVector (const DiscVecT& v, std::string varName);
    ///@}
};

/// \brief Base-class for the output of a single function in Ensight6 Case format.
///
/// We employ the command pattern: 'Describe' is the interface for registration in Ensight6OutCL.
/// 'put' is called for the output of the function at time t. The command objects are stored in Ensight6OutCL.
class Ensight6VariableCL
{
  private:
    std::string varName_,
                fileName_;
    bool        timedep_;

  public:
    Ensight6VariableCL (std::string varName, std::string fileName, bool timedep)
        : varName_( varName), fileName_( fileName), timedep_( timedep) {}
    virtual ~Ensight6VariableCL () {}

    std::string varName  () const { return varName_; }  ///< Name of the variable in einsight; also used as identifier in Ensight6OutCL.
    std::string fileName () const { return fileName_; } ///< Name of the file; for time-dependent objects, the timecode is attached by Ensight6OutCL.
    bool        Timedep  () const { return timedep_; }  ///< Is the object time-dependent?

    /// \brief Called by Ensight6OutCL::Register().
    virtual void Describe (Ensight6OutCL&) const= 0;
    /// \brief Called by Ensight6OutCL::Write().
    virtual void put      (Ensight6OutCL&) const= 0;
};

///\brief Output a geometry.
///
/// This outputs a triangulation of a multigrid.
class Ensight6GeomCL : public Ensight6VariableCL
{
  private:
    MultiGridCL* mg_;  ///< The multigrid
    int          lvl_; ///< Level of the triangulation

  public:
    Ensight6GeomCL (MultiGridCL& mg, int lvl, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), mg_( &mg), lvl_( lvl) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeGeom( this->varName());  }
    void put      (Ensight6OutCL& cf) const { cf.putGeom( *mg_, lvl_, varName()); }
};

///\brief Create an Ensight6GeomCL with operator new.
///
/// This is just for uniform code; the analoguous functions for scalars and vectors are more useful because
/// they help to avoid template parameters in user code.
inline Ensight6GeomCL&
make_Ensight6Geom (MultiGridCL& mg, int lvl, std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6GeomCL( mg, lvl, varName, fileName, timedep);
}

///\brief Represents a scalar Drops-function (P1 or P2, given as PXEvalCL) as Ensight6 variable.
template <class DiscScalarT>
class Ensight6ScalarCL : public Ensight6VariableCL
{
  private:
    const DiscScalarT f_;

  public:
    Ensight6ScalarCL (const DiscScalarT& f, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), f_( f) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeVariable( this->varName(), true); }
    void put      (Ensight6OutCL& cf) const { cf.putScalar( f_, varName()); }
};

///\brief Create an Ensight6ScalarCL<> with operator new.
///
/// This function does the template parameter deduction for user code.
template <class DiscScalarT>
  Ensight6ScalarCL<DiscScalarT>&
    make_Ensight6Scalar (const DiscScalarT& f, std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6ScalarCL<DiscScalarT>( f, varName, fileName, timedep);
}

///\brief Represents a vector Drops-function (P1 or P2, given as PXEvalCL) as Ensight6 variable.
template <class DiscVectorT>
class Ensight6VectorCL : public Ensight6VariableCL
{
  private:
    const DiscVectorT f_;

  public:
    Ensight6VectorCL (const DiscVectorT& f, std::string varName, std::string fileName, bool timedep= false)
        : Ensight6VariableCL( varName, fileName, timedep), f_( f) {}

    void Describe (Ensight6OutCL& cf) const { cf.DescribeVariable( this->varName(), false); }
    void put      (Ensight6OutCL& cf) const { cf.putVector( f_, varName()); }
};

///\brief Create an Ensight6VectorCL<> with operator new.
///
/// This function does the template parameter deduction for user code.
template <class DiscVectorT>
  Ensight6VectorCL<DiscVectorT>&
    make_Ensight6Vector (const DiscVectorT& f, std::string varName, std::string fileName, bool timedep= false)
{
    return *new Ensight6VectorCL<DiscVectorT>( f, varName, fileName, timedep);
}

#else   // parallel implementation

/// \brief Class for writing out results of a simulation in Ensight6 Case format
class EnsightP2SolOutCL
/** Specify in the constructor of this class, if the output should be in ASCII
    or binary format. In the parallel version of DROPS you can also specify if
    only the master should write out the data (this is the default behavior).
    If you do not choose this way, you have to transform the data by the
    external program "drops2ensight". */
{
  private:
    const MultiGridCL* _MG;
    const IdxDescCL*   _idx;
    char               _decDigits;
    Uint               _timestep, _numsteps;
    double             _lasttime;
    std::ostringstream _descstr,    // stores description info
                       _timestr;    // stores time info
    std::ofstream      _case;
    std::string        _geom;
    const bool         binary_;
    const int          _tag;            // tags used by this class (set to 2001) (used: tag, tag+1, tag+2)
    char               _procDigits;     // digits used for encoding procnumber
    Ulint              _nodes;          // number of exclusive vertices of this proc
    bool               _masterout;      // just the master processor writes a single file

    void AppendProccode( std::string&) const;
    void AppendTimecode( std::string&) const;
    void putTime( double time);
    void CheckFile( const std::ofstream&) const;

  public:
    // idx must refer to a numbering of the verts and edges of a certain
    // triangulation, i.e. use LevelsetP2CL::CreateNumbering before
    // constructing an EnsightSolOutCL-object.
    EnsightP2SolOutCL( const MultiGridCL& mg, const IdxDescCL* idx, bool binary=true, __UNUSED__ bool masterout=true)
      : _MG( &mg), _idx( idx), _decDigits( 0), _timestep( -1u), _numsteps( 0), _lasttime(-1), binary_( binary),
        _tag(2001), _nodes(0), _masterout(masterout)
    {
        _procDigits= 1;
        int procs=ProcCL::Size();
        while( procs>9)
        { ++_procDigits; procs/=10; }
    }
    ~EnsightP2SolOutCL() { if (_case.is_open()) CaseEnd(); }

    /// \brief Set up a case file
    /** call CaseBegin() before any other member function */
    void CaseBegin ( const char casefileName[], Uint numsteps= 0);
    /// \brief Describe geometry
    void DescribeGeom  ( const char geoName[], std::string fileName, bool timedep= false);
    /// \brief Describe scalar value finite element function
    void DescribeScalar( const char varName[], std::string fileName, bool timedep= false);
    /// \brief Describe vectorial value finite element function
    void DescribeVector( const char varName[], std::string fileName, bool timedep= false);
    /// \brief Write the multigrid (geometry) into a file
    void putGeom      ( std::string, double t= -1);
    /// \brief Write a scalar value finite element function into a file
    template<class DiscScalT>
    void putScalar    ( std::string, const DiscScalT&, double t= -1);
    /// \brief Write a vectorial value finite element function into a file
    template<class DiscVecT>
    void putVector    ( std::string, const DiscVecT&, double t= -1);
    /// \brief Rewrites case file
    void Commit    ();
    /// \brief End writing ensight files
    /** call CaseEnd() after finishing all other output */
    void CaseEnd   ();

    /// \brief Set format, that each process writes out only own values
    void SetMultipleOut() { _masterout=false; }
    /// \brief Set format, that master process writes out all values
    void SetMasterOut() { _masterout=true; }
};

//=====================================================
// Derived classes for easier usage
//=====================================================
/// \brief Class for writing out 2 phase flows in ensight format
template<typename StokesT, typename LevelsetT>
class Ensight2PhaseOutCL : public EnsightP2SolOutCL
{
  private:
    typedef EnsightP2SolOutCL base_;

  private:
    std::string file_geo_, file_scl_, file_pr_, file_vel_;  // filenames
    bool adaptive_;                                         // changing geometry
    bool timedep_;                                          // time-dependent problem
    const StokesT&   stokes_;                               // (Navier-)Stokes problem
    const LevelsetT& lset_;                                 // Levelset problem

  public:
    /**
       \brief Construct a class for writing out a 2 phase flow problem in ensight format

       All information to this class is given by this constructor:
       \param mg Multigrid of the problem
       \param idx P2 index class
       \param stokes Stokes problem class
       \param lset   Levelset problem class
       \param directory directory of ensight files
       \param caseName  name of the case file
       \param geomName name of the geometry
       \param numsteps number of time steps
       \param adaptive flag if the geometry changes over time
       \param binary    binary of ascii output
       \param masterout (in parallel) flag if the output should be done by master
    */

    Ensight2PhaseOutCL( const MultiGridCL& mg, const IdxDescCL* idx,
                        const StokesT& stokes, const LevelsetT& lset,
                        const std::string& directory, const std::string& caseName,
                        const std::string& geomName, bool adaptive,
                        Uint numsteps=0, bool binary=false, bool masterout=true )
        : base_(mg, idx, binary, masterout), adaptive_(adaptive),
          timedep_(numsteps>0), stokes_(stokes), lset_(lset)
    {
        const std::string base_filename= directory + "/" + caseName,
                          case_filename= caseName + ".case";
        file_geo_= base_filename + ".geo";
        file_scl_= base_filename + ".scl";
        file_pr_ = base_filename + ".pr";
        file_vel_= base_filename + ".vel";

        base_::CaseBegin( case_filename.c_str(), numsteps);
        base_::DescribeGeom  ( geomName.c_str(), file_geo_, adaptive_);
        base_::DescribeScalar( "Levelset",       file_scl_, timedep_);
        base_::DescribeScalar( "Pressure",       file_pr_,  timedep_);
        base_::DescribeVector( "Velocity",       file_vel_, timedep_);
    }

    ~Ensight2PhaseOutCL() { base_::CaseEnd(); }

    /// \brief Write ensight files for a time step
    void write(){
        double time   = timedep_  ? stokes_.t : -1.;
        double timeGeo= adaptive_ ? stokes_.t : -1.;
        base_::putGeom  ( file_geo_, timeGeo);
        base_::putScalar( file_scl_, lset_.GetSolution(),      time);
        base_::putScalar( file_pr_,  stokes_.GetPrSolution(),  time);
        base_::putVector( file_vel_, stokes_.GetVelSolution(), time);
        base_::Commit();
    }
};

#endif          // end of parallel implementation

class ReadEnsightP2SolCL
// read solution from Ensight6 Case format
{
  private:
    const MultiGridCL* _MG;
    const bool         binary_;
#ifdef _PAR
    char               _procDigits;     // digits used for encoding procnumber

    void AppendProccode( std::string&) const;
#endif

    void CheckFile( const std::ifstream&) const;

  public:
    ReadEnsightP2SolCL( const MultiGridCL& mg, bool binary=true)
      : _MG(&mg), binary_(binary)
    {
#ifdef _PAR
        _procDigits= 1;
        int procs=ProcCL::Size();
        while( procs>9)
        { ++_procDigits; procs/=10; }
#endif
    }

    template<class BndT>
    void ReadScalar( const std::string&, VecDescCL&, const BndT&) const;
    template<class BndT>
    void ReadVector( const std::string&, VecDescCL&, const BndT&) const;
};


//=====================================================
//              template definitions
//=====================================================

#ifndef _PAR

template<class DiscScalT>
void Ensight6OutCL::putScalar (const DiscScalT& v, std::string varName)
{
    const MultiGridCL& mg= v.GetMG();
    const Uint lvl= v.GetLevel();
    char buffer[80];
    std::memset(buffer,0,80);
    showFloat sFlo;

    std::ofstream os;
    OpenFile( os, varName);

    v.SetTime( time_);

    if(binary_)
    {
        std::strcpy(buffer,"DROPS data file, scalar variable:");
        os.write(buffer,80);

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            sFlo.f=v.val(*it);
            os.write(sFlo.s,sizeof(float));
        }

        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            sFlo.f=v.val(*it,0.5);
            os.write(sFlo.s,sizeof(float));
        }
    }
    else //ASCII-Ausgabe
    {
        int cnt=0;
        os.flags(std::ios_base::scientific);
        os.precision(5);
        os.width(12);

        os << "DROPS data file, scalar variable:\n";
        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            os << std::setw(12) << v.val( *it);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }

        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            os << std::setw(12) << v.val( *it, 0.5);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }
        os << '\n';
    }
}

template<class DiscVecT>
void Ensight6OutCL::putVector (const DiscVecT& v, std::string varName)
{
    const MultiGridCL& mg= v.GetMG();
    const Uint lvl= v.GetLevel();
    char buffer[80];
    std::memset(buffer,0,80);
    showFloat sFlo;

    std::ofstream os;
    OpenFile( os, varName);

    v.SetTime( time_);

    if(binary_)
    {
        std::strcpy(buffer,"DROPS data file, vector variable:");
        os.write( buffer, 80);

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            for (int i=0; i<3; ++i)
            {
                sFlo.f=v.val( *it)[i];
                os.write(sFlo.s,sizeof(float));
            }
        }
        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            for (int i=0; i<3; ++i)
            {
                 sFlo.f=v.val( *it)[i];
                 os.write(sFlo.s,sizeof(float));
            }
        }
    }
    else
    { // ASCII
        int cnt=0;
        os.flags(std::ios_base::scientific);
        os.precision(5);
        os.width(12);

        os << "DROPS data file, vector variable:\n";

        DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
            for (int i=0; i<3; ++i)
                os << std::setw(12) << v.val( *it)[i];
            if ( (cnt+=3)==6)
            { // Ensight expects six real numbers per line
                 cnt= 0;
                 os << '\n';
            }
        }
        DROPS_FOR_TRIANG_CONST_EDGE( mg, lvl, it) {
            for (int i=0; i<3; ++i)
                os << std::setw(12) << v.val( *it)[i];
            if ( (cnt+=3)==6)
            { // Ensight expects six real numbers per line
                 cnt= 0;
                 os << '\n';
            }
        }
        os << '\n';
    }
}

#else           // parallel implementation

template<class DiscScalT>
void EnsightP2SolOutCL::putScalar( std::string fileName, const DiscScalT& v, double t)
/** This function writes out the values of a finite element function on a
    multigrid in ensight file format. Therefore it writes out the values on
    vertices and edges in the same ordering as the function putGeom. Note, call
    DescribeScalar() and CaseBegin() before calling this function the first
    time. The template parameter DiscScalT describe the finite element function.
    \param fileName name of the file
    \param v        finite element function
    \param t        discrete time
*/
{
    const Uint lvl= _idx->TriangLevel();
    char buffer[80];
    showFloat sFlo;

    if ( t!=-1){
        putTime( t);
        AppendTimecode( fileName);
    }
    if (!_masterout)
        AppendProccode(fileName);

    std::ofstream os( fileName.c_str());
    CheckFile( os);

    // parallele Ausgabe
    int cnt=0;
    if(!binary_)
    {
        os.flags(std::ios_base::scientific);
        os.precision(5);
        os.width(12);
    }
    VectorCL *vals=0;

    if (!_masterout)
    {
        os << _nodes << "  ";
        os << "DROPS data file, scalar variable:\n";
    }
    else
    {
        vals = new VectorCL(_nodes);
        if (ProcCL::IamMaster()) {
            if(binary_) {
                std::strcpy(buffer,"DROPS data file, scalar variable:");
                os.write(buffer,80);
            }
            else {
                os << "DROPS data file, scalar variable:\n";
            }
        }
    }

    Uint pos=0;

    for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
        end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        // Write only the vals on exclusive vertices!
        if (!it->IsExclusive(PrioHasUnk))
            continue;

        if (!_masterout)
        {
            os << std::setw(12) << v.val( *it);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }
        else
        {
            (*vals)[pos++]=v.val( *it);
        }
    }

    for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        // Write only the vals on exclusive edges!
        if (!it->IsExclusive(PrioHasUnk))
            continue;

        if (!_masterout)
        {
            os << std::setw(12) << v.val( *it, 0.5);
            if ( (++cnt)==6)
            { // Ensight expects six real numbers per line
                cnt= 0;
                os << '\n';
            }
        }
        else
        {
            (*vals)[pos++]=v.val( *it, 0.5);
        }
    }

    if (_masterout)
    {
        if(pos!=vals->size())
            throw DROPSErrCL("EnsightP2SolOutCL::putScalar: Not enough data collected!");

        Uint maxUnk= ProcCL::GlobalMax(_nodes,ProcCL::Master());

        if (!ProcCL::IamMaster())
        {
            ProcCL::RequestT req;
            req= ProcCL::Isend(Addr(*vals), _nodes, ProcCL::Master(), _tag);
            ProcCL::Wait(req);
        }
        else
        {
            VectorCL *recvBuf = new VectorCL(maxUnk);
            VectorCL *tmpBuf=0;
            IdxT numUnk;

            for (int p=0; p<ProcCL::Size(); ++p)
            {
                if (p!=ProcCL::MyRank())
                {
                    ProcCL::StatusT stat;
                    ProcCL::Probe(p, _tag, stat);
                    numUnk = (IdxT)ProcCL::GetCount<double>(stat);

                    ProcCL::Recv(Addr(*recvBuf), numUnk, p, _tag);
                }
                else
                {
                    tmpBuf=recvBuf;
                    recvBuf=vals;
                    numUnk=_nodes;
                }
                for (IdxT i=0; i<numUnk; ++i)
                {
                    if(binary_) {
                        sFlo.f=(*recvBuf)[i];
                        os.write(sFlo.s,sizeof(float));
                    }
                    else {
                        os << std::setw(12) << (*recvBuf)[i];
                        if ( (++cnt)==6)
                        { // Ensight expects six real numbers per line
                            cnt= 0;
                            os << '\n';
                        }
                    }
                }
                if (p==ProcCL::MyRank())
                {
                    recvBuf=tmpBuf;
                    tmpBuf=0;
                }
            }
            delete recvBuf;
            if(!binary_) {
                os << '\n';
            }
        }
        delete vals;
    }
    else {
        if(!binary_) {
            os << '\n';
        }
    }
    os.close();
}

template<class DiscVecT>
void EnsightP2SolOutCL::putVector( std::string fileName, const DiscVecT& v, double t)
/** This function writes out the values of a finite element function on a
    multigrid in ensight file format. Therefore it writes out the values on
    vertices and edges in the same ordering as the function putGeom. Note, call
    DescribeScalar() and CaseBegin() before calling this function the first
    time. The template parameter DiscScalT describe the finite element function.
    \param fileName name of the file
    \param v        finite element function
    \param t        discrete time
*/
{
    const Uint lvl= _idx->TriangLevel();

    if ( t!=-1){
        putTime( t);
        AppendTimecode( fileName);
    }
    if (!_masterout)
        AppendProccode(fileName);

    std::ofstream os( fileName.c_str());
    CheckFile( os);
    //os.flags(std::ios_base::scientific);

    if(binary_)
    {
        char buffer[128];
        showFloat sFlo;

        // parallele Ausgabe
        VectorCL *vals=0;

        if (_masterout)
        {
            vals = new VectorCL(3*_nodes);
            if (ProcCL::IamMaster())
            {
                std::strcpy(buffer,"DROPS data file, vector variable:");
                os.write(buffer,80);
            }
        }

        Uint pos=0;

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {

            // Write only the vals on exclusive vertices!
            if (!it->IsExclusive(PrioHasUnk))
                continue;
            if (_masterout)
            {
                for (int i=0; i<3; ++i)
                    (*vals)[pos++]= v.val( *it)[i];
            }
        }
        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            // Write only the vals on exclusive edges!
            if (!it->IsExclusive(PrioHasUnk))
                continue;

            if (_masterout)
            {
                for (int i=0; i<3; ++i)
                    (*vals)[pos++]= v.val( *it)[i];
            }
        }

        if (_masterout)
        {
            if(pos!=vals->size())
                throw DROPSErrCL("EnsightP2SolOutCL::putVector: Not enough data collected!");

            Uint maxUnk= ProcCL::GlobalMax(3*_nodes, ProcCL::Master());

            if (!ProcCL::IamMaster())
            {
                ProcCL::RequestT req;
                req= ProcCL::Isend(Addr(*vals), 3*_nodes, ProcCL::Master(), _tag);
                ProcCL::Wait(req);
            }
            else
            {
                VectorCL *recvBuf = new VectorCL(maxUnk);
                VectorCL *tmpBuf=0;
                IdxT numUnk;

                for (int p=0; p<ProcCL::Size(); ++p)
                {
                    if (p!=ProcCL::MyRank())
                    {
                        ProcCL::StatusT stat;
                        ProcCL::Probe(p, _tag, stat);
                        numUnk = (IdxT)ProcCL::GetCount<double>(stat);

                        ProcCL::Recv(Addr(*recvBuf), numUnk, p, _tag);
                    }
                    else
                    {
                        tmpBuf=recvBuf;
                        recvBuf=vals;
                        numUnk=3*_nodes;
                    }

                    for (IdxT i=0; i<numUnk; ++i)
                    {
                     //   os << std::setw(12) << (*recvBuf)[i];
                        sFlo.f =(float) (*recvBuf)[i];
                        os.write(sFlo.s,sizeof(float));

                    }
                    if (p==ProcCL::MyRank())
                    {
                        recvBuf=tmpBuf;
                        tmpBuf=0;
                    }
                }
                delete recvBuf;

            }
            delete vals;
        }

    }
    else
    { // ASCII
        int cnt=0;
        os.precision(5);
        os.width(12);

        VectorCL *vals=0;

        if (!_masterout)
        {
            os << (3*_nodes) << "  ";
            os << "DROPS data file, vector variable:\n";
        }
        else
        {
            vals = new VectorCL(3*_nodes);
            if (ProcCL::IamMaster())
                os << "DROPS data file, vector variable:\n";
        }

        Uint pos=0;

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {

            // Write only the vals on exclusive vertices!
            if (!it->IsExclusive(PrioHasUnk))
                continue;
            if (!_masterout)
            {
                for (int i=0; i<3; ++i)
                    os << std::setw(12) << v.val( *it)[i];
                if ( (cnt+=3)==6)
                { // Ensight expects six real numbers per line
                    cnt= 0;
                    os << '\n';
                }
            }
            else
            {
                for (int i=0; i<3; ++i)
                    (*vals)[pos++]= v.val( *it)[i];
            }
        }
        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            // Write only the vals on exclusive edges!
            if (!it->IsExclusive(PrioHasUnk))
                continue;

            if (!_masterout)
            {
                for (int i=0; i<3; ++i)
                    os << std::setw(12) << v.val( *it)[i];
                if ( (cnt+=3)==6)
                { // Ensight expects six real numbers per line
                    cnt= 0;
                    os << '\n';
                }
            }
            else
            {
                for (int i=0; i<3; ++i)
                    (*vals)[pos++]= v.val( *it)[i];
            }
        }

        if (_masterout)
        {
            if(pos!=vals->size())
                throw DROPSErrCL("EnsightP2SolOutCL::putVector: Not enough data collected!");

            Uint maxUnk= ProcCL::GlobalMax(3*_nodes, ProcCL::Master());

            if (!ProcCL::IamMaster())
            {
                ProcCL::RequestT req;
                req= ProcCL::Isend(Addr(*vals), 3*_nodes, ProcCL::Master(), _tag);
                ProcCL::Wait(req);
            }
            else
            {
                VectorCL *recvBuf = new VectorCL(maxUnk);
                VectorCL *tmpBuf=0;
                IdxT numUnk;

                for (int p=0; p<ProcCL::Size(); ++p)
                {
                    if (p!=ProcCL::MyRank())
                    {
                        ProcCL::StatusT stat;
                        ProcCL::Probe(p, _tag, stat);
                        numUnk = (IdxT)ProcCL::GetCount<double>(stat);

                        ProcCL::Recv(Addr(*recvBuf), numUnk, p, _tag);
                    }
                    else
                    {
                        tmpBuf=recvBuf;
                        recvBuf=vals;
                        numUnk=3*_nodes;
                    }

                    for (IdxT i=0; i<numUnk; ++i)
                    {
                        os << std::setw(12) << (*recvBuf)[i];
                        if ( (++cnt)==6)
                        { // Ensight expects six real numbers per line
                            cnt= 0;
                            os << '\n';
                        }
                    }
                    if (p==ProcCL::MyRank())
                    {
                        recvBuf=tmpBuf;
                        tmpBuf=0;
                    }
                }
                delete recvBuf;
                os << '\n';
            }
            delete vals;
        }
        else
            os << '\n';
    }
    os.close();
}

#endif          // end of parallel implementation


// ========== ReadEnsightP2SolCL ==========

template <class BndT>
void ReadEnsightP2SolCL::ReadScalar( const std::string& file, VecDescCL& v, const BndT& bnd) const
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
    std::string fileName(file);
#ifdef _PAR
    const ExchangeCL& ex= v.RowIdx->GetEx();
    AppendProccode(fileName);
#endif

    std::ifstream is( fileName.c_str());
    CheckFile( is);

    if (binary_)
    {
#ifdef _PAR
        IF_MASTER
            std::cout << "ReadEnsightP2SolCL::ReadScalar: Binary parallel output for ensight not yet implemented! (of)" << std::endl;
        throw DROPSErrCL("ReadEnsightP2SolCL::ReadScalar: Binary parallel output for ensight not yet implemented!");
#endif
        showFloat fl;
        char buffer[80];
        is.read( buffer, 80);           //ignore first 80 characters

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= (double)fl.f;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= fl.f;
        }
    }
    else
    { // ASCII
        char buf[256];
        double d= 0;

        is.getline( buf, 256); // ignore first line

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
             end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
#ifdef _PAR
            if (!it->IsExclusive(PrioHasUnk)) continue;
#endif

            is >> d;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= d;
        }
        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
#ifdef _PAR
            if (!it->IsExclusive(PrioHasUnk)) continue;
#endif

            is >> d;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            v.Data[it->Unknowns(idx)]= d;
        }
#ifdef _PAR
        ex.Accumulate(v.Data);
#endif
    }

    CheckFile( is);
}

template <class BndT>
void ReadEnsightP2SolCL::ReadVector( const std::string& file, VecDescCL& v, const BndT& bnd) const
{
    const Uint lvl= v.GetLevel(),
               idx= v.RowIdx->GetIdx();
    std::string fileName(file);
#ifdef _PAR
    const ExchangeCL& ex= v.RowIdx->GetEx();
    AppendProccode(fileName);
#endif

    std::ifstream is( fileName.c_str());
    CheckFile( is);

    double d0= 0, d1= 0, d2= 0;

    if(binary_)
    {
#ifdef _PAR
        IF_MASTER
            std::cout << "ReadEnsightP2SolCL::ReadVector: Binary parallel output for ensight not yet implemented! (of)" << std::endl;
        throw DROPSErrCL("ReadEnsightP2SolCL::ReadVector: Binary parallel output for ensight not yet implemented!");
#endif
        showFloat fl;
        //std::cout<<"READVECTOR: "<<file.c_str()<<"\n";
        char buffer[80];
        is.read( buffer, 80);       //ignore first 80 characters

        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            d0=fl.f;
            is.read( fl.s, sizeof(float));
            d1=fl.f;
            is.read( fl.s, sizeof(float));
            d2=fl.f;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
             end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
            is.read( fl.s, sizeof(float));
            d0=fl.f;
            is.read( fl.s, sizeof(float));
            d1=fl.f;
            is.read( fl.s, sizeof(float));
            d2=fl.f;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }
    }
    else
    { // ASCII
        char buf[256];

        is.getline( buf, 256); // ignore first line

        int count=0;
        for (MultiGridCL::const_TriangVertexIteratorCL it= _MG->GetTriangVertexBegin(lvl),
            end= _MG->GetTriangVertexEnd(lvl); it!=end; ++it)
        {
#ifdef _PAR
            if (!it->IsExclusive(PrioHasUnk)) continue;
#endif
            is >> d0 >> d1 >> d2;
            count +=3;

            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;
            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }

        for (MultiGridCL::const_TriangEdgeIteratorCL it= _MG->GetTriangEdgeBegin(lvl),
            end= _MG->GetTriangEdgeEnd(lvl); it!=end; ++it)
        {
#ifdef _PAR
            if (!it->IsExclusive(PrioHasUnk)) continue;
#endif
            is >> d0 >> d1 >> d2;
            count +=3;
            if (bnd.IsOnDirBnd( *it) || !(it->Unknowns.Exist(idx)) ) continue;

            const IdxT Nr= it->Unknowns(idx);
            v.Data[Nr]= d0;    v.Data[Nr+1]= d1;    v.Data[Nr+2]= d2;
        }
#ifdef _PAR
        ex.Accumulate(v.Data);

        if (!is)
            std::cout << "["<<ProcCL::MyRank()<<"] Read "<<count<<" numbers from file <"<<fileName<<">"<<std::endl;
#endif
    }

    CheckFile( is);
}

} // end of namespace DROPS

#endif
