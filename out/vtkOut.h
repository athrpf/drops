/// \file vtkOut.h
/// \brief solution output in VTK legacy format
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

class VTKVariableCL; //forward declaration

/// \brief Class for writing out results of a simulation in VTK legacy format
class VTKOutCL
/** This class writes out data in VTK legacy format. Therefore all data are sent
    to the master processor, that is responsible for the I/O. The user can write
    out the geometry, scalar and vector-valued finite element functions.

    To use this class correct, the geometry has to be written out first. Then
    the finite element functions can be put into the vtk files.

    \todo binary output does not work correctly
*/
{
  private:
#ifndef _PAR
    typedef std::map<const VertexCL*, Uint> vertexAddressMapT;
    typedef std::map<const EdgeCL*, Uint>     edgeAddressMapT;
    typedef VectorBaseCL<Uint>              TetraVecT;
#else
    typedef std::map < GIDT, Uint > GIDMapT;     // map GID to number
    typedef std::vector< Uint >        ProcOffsetCT;// position in tetras_ where all tetras of a processor are lying
    typedef VectorBaseCL<GIDT>      TetraVecT;
#endif

    const MultiGridCL& mg_;                         // reference to the multigrid
    char               decDigits_;                  // number of digits for encoding time in filename
    Uint               timestep_, numsteps_;        // actual timestep and number of timesteps
    std::string        descstr_;                    // stores description info
    std::string        filename_;                   // filenames
    std::ofstream      file_;                       // actual file where to put data
    std::map<std::string, VTKVariableCL*> vars_;    ///< The variables stored by varName.
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
    TetraVecT             tetras_;                  // Connectivities (tetras)
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
    void GatherCoord( VectorBaseCL<GIDT>&, VectorBaseCL<float>&);
    /// Communicate coords to master
    void CommunicateCoords(const VectorBaseCL<GIDT>&, const VectorBaseCL<float>&);
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
    void GatherTetra( VectorBaseCL<GIDT>&) const;
    /// Communicate tetras
    void CommunicateTetra(const VectorBaseCL<GIDT>&);
    /// Write out the distribution of tetrahedra as CELL_DATA
    void WriteDistribution();
#endif
    /// Write Tetras
    void WriteTetra();
    //@}

    /// \name Write out FE functions
    //@{
    /// Gather scalar data
    template <typename DiscScalT>
    void GatherScalar(const DiscScalT&, VectorBaseCL<float>&) const;
    /// Gather vectorial data
    template <typename DiscVecT>
    void GatherVector(const DiscVecT&, VectorBaseCL<float>&) const;
    /// Write data
    void WriteValues(const VectorBaseCL<float>&, const std::string&, int);
#ifdef _PAR
    /// Communicate data
    void CommunicateValues(const VectorBaseCL<float>&, VectorBaseCL<float>&, int);
#endif
    //@}

  public:
    /// \brief Constructor of this class
    VTKOutCL(const MultiGridCL& mg, const std::string& dataname, Uint numsteps,
             const std::string& filename, bool binary);

    /// \brief Register a variable or the geometry for output with Write().
    ///
    /// The class takes ownership of the objects, i. e. it destroys them with delete in its destructor.
    void Register (VTKVariableCL& var);
    /// Write out all registered objects.
    void Write (double time, bool writeDistribution=false);
    /// \brief Put geometry into the vtk file
    void PutGeom( double time, bool writeDistribution=false);

    /// \brief Write scalar variable into file
    template <typename DiscScalT>
    void PutScalar( const DiscScalT&, const std::string&);

    /// \brief Write vector-valued variable into file
    template <typename DiscVecT>
    void PutVector( const DiscVecT&, const std::string&);

    /// \brief End output of a file
    void Commit(){
        file_.close();
        timestep_++;
    }

    /// \brief Clear all internal data
    void Clear();
};

/// \brief Base-class for the output of a single function in Ensight6 Case format.
///
/// 'put' is called for the output of the function at time t. The command objects are stored in VTKOutCL.
class VTKVariableCL
{
  private:
    std::string varName_;

  public:
    VTKVariableCL (std::string varName)
        : varName_( varName) {}
    virtual ~VTKVariableCL () {}

    std::string varName() const { return varName_; }  ///< Name of the variable in VTK; also used as identifier in VTKOutCL.

    /// \brief Called by VTKOutCL::Write().
    virtual void put( VTKOutCL&) const= 0;
};

///\brief Represents a scalar Drops-function (P1 or P2, given as PXEvalCL) as VTK variable.
template <class DiscScalarT>
class VTKScalarCL : public VTKVariableCL
{
  private:
    const DiscScalarT f_;

  public:
    VTKScalarCL( const DiscScalarT& f, std::string varName)
        : VTKVariableCL( varName), f_( f) {}

    void put( VTKOutCL& cf) const { cf.PutScalar( f_, varName()); }
};

///\brief Create an VTKScalarCL<> with operator new.
///
/// This function does the template parameter deduction for user code.
template <class DiscScalarT>
  VTKScalarCL<DiscScalarT>&
    make_VTKScalar( const DiscScalarT& f, std::string varName)
{
    return *new VTKScalarCL<DiscScalarT>( f, varName);
}

///\brief Represents a vector Drops-function (P1 or P2, given as PXEvalCL) as VTK variable.
template <class DiscVectorT>
class VTKVectorCL : public VTKVariableCL
{
  private:
    const DiscVectorT f_;

  public:
    VTKVectorCL( const DiscVectorT& f, std::string varName)
        : VTKVariableCL( varName), f_( f) {}

    void put( VTKOutCL& cf) const { cf.PutVector( f_, varName()); }
};

///\brief Create an VTKVectorCL<> with operator new.
///
/// This function does the template parameter deduction for user code.
template <class DiscVectorT>
  VTKVectorCL<DiscVectorT>&
    make_VTKVector( const DiscVectorT& f, std::string varName)
{
    return *new VTKVectorCL<DiscVectorT>( f, varName);
}

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

    /// \brief Write all information of a two phase problem
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

template <typename DiscScalT>
  void VTKOutCL::PutScalar( const DiscScalT& f, const std::string& name)
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

template <typename  DiscVecT>
  void VTKOutCL::PutVector( const DiscVecT& f, const std::string& name)
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

template <typename DiscScalT>
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

template <typename DiscVecT>
  void VTKOutCL::GatherVector(const DiscVecT& f, VectorBaseCL<float>& locData) const
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

} // end of namespace DROPS

#endif
