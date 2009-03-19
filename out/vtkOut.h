//**************************************************************************
// File:    vtkOut.h                                                       *
// Content: solution output in VTK legacy format                           *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
//**************************************************************************

#ifndef DROPS_VTKOUT_H
#define DROPS_VTKOUT_H

#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "out/ensightOut.h"
#include <map>
#include <vector>

#ifdef _PAR
# include "parallel/parallel.h"
# include "parallel/exchange.h"
#endif


namespace DROPS
{

/// \brief Class for writing out results of a simulation in VTK legacy format
class VTKOutCL
/** This class writes out data in VTK legacy format. Therefore all data are sent
    to the master processor, that is responsible for the I/O. The user can write
    out the geometry, scalar and vectorial valued finite element functions.

    To use this class correct, the geometry has to be written out at first. Then
    the finite element functions can be put into the vtk files.

    \todo binary output does not correctly work
*/
{
  private:
#ifndef _PAR
    typedef std::map<const VertexCL*, Uint> vertexAddressMapT;
    typedef std::map<const EdgeCL*, Uint>     edgeAddressMapT;
#else
    typedef std::map < DDD_GID, Uint > GIDMapT;     // map GID to number
    typedef std::vector< Uint >        ProcOffsetCT;// position in tetras_ where all tetras of a processor are lying
#endif

    const MultiGridCL& mg_;                         // reference to the multigrid
    char               decDigits_;                  // number of digits for encoding time in filename
    Uint               timestep_, numsteps_;        // actual timestep and number of timesteps
    std::string        descstr_;                    // stores description info
    std::string        filename_;                   // filenames
    std::ofstream      file_;                       // actual file where to put data
    const bool         binary_;                     // output in binary or ascii format
    bool               geomwritten_;                // flag if geometry has been written

#ifndef _PAR
    vertexAddressMapT   vAddrMap_;                  // Map vertex adress to a unique (consecutive) number
    edgeAddressMapT     eAddrMap_;                  // Map edge Adress to a unique (consecutive) number
#else
    size_t              mgVersion_;                 // version number of last written out multigrid (sequential version of DROPS do not provide a version-number for the MG)
    const int           tag_;                       // tags used by this class (set to 3001)
    GIDMapT             GIDAddrMap_;                // Map GID to unique (consecutive) number (index to a number)
    ProcOffsetCT        procOffset_;                // index in tetras_, where tetras of a processor are lying
#endif

    VectorBaseCL<float>   coords_;                  // Coordinates of the points
    VectorBaseCL<DDD_GID> tetras_;                  // Connectivities (tetras)
    Uint                  numPoints_;               // number of points (only accessible by master process)
    Uint                  numTetras_;               // number of tetras (only accessible by master process)
    Uint                  numLocPoints_;            // number of local exclusive verts and edges
    bool                  wrotePointDataLine_;      // flag if descibtionline for point data has been written

    /// Put timecode as postfix to a filename
    void AppendTimecode( std::string&) const;
    /// Check if the file is open
    void CheckFile( const std::ofstream&) const;
    /// Create new file
    void NewFile( double time);
    /// Put describtion header into file
    void PutHeader( double time);

#ifdef _PAR
    /// Check if the calling processor should gather information on the given simplex
    template<typename SimplexT>
    bool AmIResponsible(const SimplexT& s) const{
        return s.IsExclusive(PrioHasUnk);
    }
#endif

    /// \name Write out coordinates of vertices
    //@{
#ifndef _PAR
    /// Gather coordinates and (g)ids
    void GatherCoord();
#else
    /// Gather coordinates and (g)ids
    void GatherCoord( VectorBaseCL<DDD_GID>&, VectorBaseCL<float>&);
    /// Communicate coords to master
    void CommunicateCoords(const VectorBaseCL<DDD_GID>&, const VectorBaseCL<float>&);
#endif
    /// Write coordinates into a vtk legacy file
    void WriteCoords();
    //@}

    /// \name Write out connectivities
    //@{
#ifndef _PAR
    /// Gather tetras
    void GatherTetra();
#else
    /// Gather tetras
    void GatherTetra( VectorBaseCL<DDD_GID>&) const;
    /// Communicate tetras
    void CommunicateTetra(const VectorBaseCL<DDD_GID>&);
    /// Write out the distribution of tetrahedra as CELL_DATA
    void WriteDistribution();
#endif
    /// Write Tetras
    void WriteTetra();
    //@}

    /// \name Write out FE functions
    //@{
    /// Gather scalar data
    template<class DiscScalT>
    void GatherScalar(const DiscScalT&, VectorBaseCL<float>&) const;
    /// Gather vectorial data
    template<class DiscScalT>
    void GatherVector(const DiscScalT&, VectorBaseCL<float>&) const;
    /// Write data
    void WriteValues(const VectorBaseCL<float>&, const std::string&, int);
#ifdef _PAR
    /// Communicar data
    void CommunicateValues(const VectorBaseCL<float>&, VectorBaseCL<float>&, int);
#endif
    //@}

  public:
    /// \brief Constructor of this class
    VTKOutCL(const MultiGridCL& mg, const std::string& dataname, Uint numsteps,
             const std::string& filename, bool binary);

    /// \brief Put geometry into the vtk file
    void PutGeom( double time, bool writeDistribution=false);

    /// \brief Write scalar variable into file
    template<class DiscScalT>
    void PutScalar( const std::string&, const DiscScalT&);

    /// \brief Write vectorial variable into file
    template<class DiscScalT>
    void PutVector( const std::string&, const DiscScalT&);

    /// \brief End output of a file
    void Commit(){
        file_.close();
        timestep_++;
    }

    /// \brief Clear all internal data
    void Clear();
};

//=====================================================
// Derived classes for easier usage
//=====================================================

/// \brief Class for writing out a twophase flow in VTK file format
template<typename StokesCL, typename LevelsetCL>
class TwoPhaseVTKCL : public VTKOutCL
/** This class writes out information about the geometry, velocity, pressure
    and the level-set function.*/
{
  private:
    typedef VTKOutCL base_;

  protected:
    const StokesCL&   stokes_;
    const LevelsetCL& levelset_;
    const bool        writeDist_;

  public:
    /// \brief Constructor of this class
    TwoPhaseVTKCL (const MultiGridCL& mg, const StokesCL& st, const LevelsetCL& ls,
                   Uint numsteps, const std::string& filename, bool binary, bool writeDistribution=false);

    /// \brief Write alle information of a two phase problem
    void write();
};

/// \brief Class for writing out a twophase flow with transport in VTK file format
template<typename StokesCL, typename LevelsetCL, typename TransportCL>
class TwoPhaseTransportVTKCL : public TwoPhaseVTKCL<StokesCL, LevelsetCL>
/** This class writes out information about the geometry, velocity, pressure,
    the level-set function and transport .*/
{
  private:
    typedef TwoPhaseVTKCL<StokesCL, LevelsetCL> base_;

    const TransportCL& transport_;

  public:
    /// \brief Constructor of this class
    TwoPhaseTransportVTKCL (const MultiGridCL& mg, const StokesCL& st, const LevelsetCL& ls,
                            const TransportCL& tr, Uint numsteps, const std::string& filename,
                            bool binary);

    /// \brief Write alle information of a twophase problem with mass transport
    void write();
};

//=====================================================
//              template definitions
//=====================================================

template<class DiscScalT>
  void VTKOutCL::PutScalar( const std::string& name, const DiscScalT& f)
/** Write values of a scalar valued function into the vtk legacy file
    \param name name of the function
    \param f    function*/
{
    VectorBaseCL<float> allData;

#ifndef _PAR
    GatherScalar( f, allData);
#else
    VectorBaseCL<float> locData;
    GatherScalar( f, locData);
    CommunicateValues( locData, allData, 1);
#endif
    WriteValues( allData, name, 1);
}

template<class DiscScalT>
  void VTKOutCL::PutVector( const std::string& name, const DiscScalT& f)
/** Write values of a scalar valued function into the vtk legacy file
    \param name name of the function
    \param f    function*/
{
    VectorBaseCL<float> allData;

#ifndef _PAR
    GatherVector( f, allData);
#else
    VectorBaseCL<float> locData;
    GatherVector( f, locData);
    CommunicateValues( locData, allData, 3);
#endif
    WriteValues( allData, name, 3);
}

template<class DiscScalT>
  void VTKOutCL::GatherScalar(const DiscScalT& f, VectorBaseCL<float>& locData) const
/** Gather values of the function f on vertices and edges*/
{
    // Allocate mem
    locData.resize(numLocPoints_);

    // Get values of vertices
    Uint pos= 0;
    for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(); it!=mg_.GetTriangVertexEnd(); ++it){
#ifdef _PAR
        if (AmIResponsible(*it))
#endif
            locData[pos++]= (float)f.val( *it);
    }

    // Get values on edges
    for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(); it!=mg_.GetTriangEdgeEnd(); ++it){
#ifdef _PAR
        if (AmIResponsible(*it))
#endif
            locData[pos++]= (float)f.val( *it, 0.5);
    }
}

template<class DiscScalT>
  void VTKOutCL::GatherVector(const DiscScalT& f, VectorBaseCL<float>& locData) const
/** Gather values of the function f on vertices and edges*/
{
    // Allocate mem
    locData.resize(3*numLocPoints_);

    // Get values of vertices
    Uint pos= 0;
    for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(); it!=mg_.GetTriangVertexEnd(); ++it){
#ifdef _PAR
        if (AmIResponsible(*it))
#endif
            for (int j=0; j<3; ++j)
                locData[pos++]= (float)f.val( *it)[j];
    }

    // Get values on edges
    for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(); it!=mg_.GetTriangEdgeEnd(); ++it){
#ifdef _PAR
        if (AmIResponsible(*it))
#endif
            for (int j=0; j<3; ++j)
                locData[pos++]= (float)f.val( *it, 0.5)[j];
    }
}

template<typename StokesCL, typename LevelsetCL>
TwoPhaseVTKCL<StokesCL, LevelsetCL>::TwoPhaseVTKCL (const MultiGridCL& mg, const StokesCL& st, const LevelsetCL& ls,
  Uint numsteps, const std::string& filename, bool binary, bool writeDistribution)
/** This constructor init the base class VTKOutCL, in order to write out
    twophase flows
\param mg geometry of the given problem
\param st Stokes (or Nabier-Stokes) problem class
\param ls levelset function
\param numsteps number of timesteps, that should be written out
\param filename prefix of all files
\param binary    Write out files in binary format. (Does not work)
\param writeDistribution Write out distribution
\todo Support of binary output
*/
    : base_(mg, "DROPS data", numsteps, filename, binary),
        stokes_(st), levelset_(ls), writeDist_(writeDistribution) {}


template<typename StokesCL, typename LevelsetCL>
  void TwoPhaseVTKCL<StokesCL, LevelsetCL>::write()
/** Writes out the geometry information as well as the pressure and velocity of
    the stokes problem and the level-set function in a VTK legacy format.
*/
{
    base_::PutGeom(stokes_.t, writeDist_);
    base_::PutVector( "velocity",  stokes_.GetVelSolution() );
    base_::PutScalar( "pressure",  stokes_.GetPrSolution() );
    base_::PutScalar( "level-set", levelset_.GetSolution() );
    base_::Commit();
}


template<typename StokesCL, typename LevelsetCL, typename TransportCL>
TwoPhaseTransportVTKCL<StokesCL, LevelsetCL, TransportCL>::TwoPhaseTransportVTKCL(
    const MultiGridCL& mg, const StokesCL& st, const LevelsetCL& ls, const TransportCL& tr,
    Uint numsteps, const std::string& filename, bool binary)
/** This consstructor init the base class TwoPhaseVTKCL in order to write out
    a twophase flow with mass transport.
\param mg geometry of the given problem
\param st Stokes (or Nabier-Stokes) problem class
\param ls levelset function
\param tr transport
\param numsteps number of timesteps, that should be written out
\param filename prefix of all files
\param binary    Write out files in binary format. (Does not work)
\todo Support of binary output
*/
    : base_(mg, st, ls, numsteps, filename, binary), transport_(tr) {}


template<typename StokesCL, typename LevelsetCL, typename TransportCL>
  void TwoPhaseTransportVTKCL<StokesCL, LevelsetCL, TransportCL>::write()
/** Writes out the geometry information as well as the pressure and velocity of
    the stokes problem, the level-set function and the transport solution in a
    VTK legacy format.
*/
{
    base_::PutGeom(base_::stokes_.t);
    base_::PutVector( "velocity",  base_::stokes_.GetVelSolution() );
    base_::PutScalar( "pressure",  base_::stokes_.GetPrSolution() );
    base_::PutScalar( "level-set", base_::levelset_.GetSolution() );
    base_::PutScalar( "transport", transport_.GetSolution() );
    base_::Commit();
}

} // end of namespace DROPS

#endif
