/// \file vtkOut.h
/// \brief solution output in VTK XML format
/// \author LNM RWTH Aachen: Jens Berger, Patrick Esser, Joerg Grande, Sven Gross, Martin Horsky, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

/// \brief Class for writing out results of a simulation in VTK XML format
class VTKOutCL
/** This class writes out data in VTK XML format. The user can write
    out the geometry, scalar and vector-valued finite element functions.
*/
{
  private:

    typedef std::map<const VertexCL*, Uint>     vertexAddressMapT;
    typedef std::map<const EdgeCL*, Uint>         edgeAddressMapT;
    typedef std::map<std::string, VTKVariableCL*> VTKvarMapT;
    typedef VectorBaseCL<Uint>                    TetraVecT;

    const MultiGridCL& mg_;                         ///< reference to the multigrid
    char               decDigits_;                  ///< number of digits for encoding time in filename
    Uint               timestep_, numsteps_;        ///< actual timestep and number of timesteps
    std::string        descstr_;                    ///< stores description info
    std::string        filename_;                   ///< filenames
    std::ofstream      file_;                       ///< actual file where to put data
    VTKvarMapT         vars_;                       ///< The variables stored by varName.
    const bool         binary_;                     ///< output in binary or ascii format
    bool               geomwritten_;                ///< flag if geometry has been written

    vertexAddressMapT   vAddrMap_;                  ///< Map vertex address to a unique (consecutive) number
    edgeAddressMapT     eAddrMap_;                  ///< Map edge address to a unique (consecutive) number

    VectorBaseCL<float>   coords_;                  ///< Coordinates of the points
    TetraVecT             tetras_;                  ///< Connectivities (tetras)
    Uint                  numPoints_;               ///< number of points (only accessible by master process)
    Uint                  numTetras_;               ///< number of tetras (only accessible by master process)
    Uint                  numLocPoints_;            ///< number of local exclusive verts and edges
    bool                  wrotePointDataLine_;      ///< flag if description line for point data has been written

    /// Puts time-code as a post-fix to the filename
    void AppendTimecode( std::string&) const;
    /// Checks whether the file is open
    void CheckFile( const std::ofstream&) const;
    /// Creates new file
    void NewFile( /*TODO double time */bool writeDistribution);
    /// Puts the description header into the file
    void PutHeader();
    /// Puts the footer into the file
    void PutFooter( );
    /// Writes the variable names of the numerical data into the file
    void WriteVarNames(std::ofstream&,bool masterfile=0);

       /// \name Writes out coordinates of vertices
    //@{
    /// Gathers coordinates
    void GatherCoord();
    /// Writes coordinates into a VTK file
    void WriteCoords();
    //@}

    /// \name Write out connectivities
    //@{
    /// Gather tetras
    void GatherTetra();

    /// Writes out the distribution of the tetrahedra as CellData
    void WriteDistribution();

    /// Write Tetras
    void WriteTetra( bool writeDistribution);
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
    void WriteValues(const VectorBaseCL<float>&, const std::string&, int, std::ofstream* masterFile= 0);
    //@}



  public:
    /// \brief Constructor of this class
    VTKOutCL(const MultiGridCL& mg, const std::string& dataname, Uint numsteps,
             const std::string& filename, bool binary);
    ~VTKOutCL();

    /// \brief Register a variable or the geometry for output with Write().
    ///
    /// The class takes ownership of the objects, i. e. it destroys them with delete in its destructor.
    void Register ( VTKVariableCL& var);
    /// Writes out all registered objects.
    void Write ( double time, bool writeDistribution=true);
    /// \brief Puts geometry into the VTK file
    void PutGeom( bool writeDistribution=true);

    /// \brief Writes scalar-valued variable into the file
    template <typename DiscScalT>
    void PutScalar( const DiscScalT&, const std::string&);

    /// \brief Writes vector-valued variable into the file
    template <typename DiscVecT>
    void PutVector( const DiscVecT&, const std::string&);

    /// \brief Ends output of a file
    void Commit(){
    	PutFooter();
        file_.close();
        timestep_++;
    }

    /// \brief Clear all internal data
    void Clear();
};

/// \brief Base-class for the output of a single function in VTK format.
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
    virtual Uint GetDim() const = 0;
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
    Uint GetDim() const { return 1; }
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
    Uint GetDim() const { return 3; }
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

    /// \brief Write all information of a two-phase problem with mass transport
    void write();
};

//=====================================================
//              template definitions
//=====================================================

template <typename DiscScalT>
  void VTKOutCL::PutScalar( const DiscScalT& f, const std::string& name)
/** Writes the values of a scalar valued function into the VTK file
    \param name name of the function
    \param f    function*/
{
    VectorBaseCL<float> allData;
    GatherScalar( f, allData);
    WriteValues(allData, name, 1);
}

template <typename  DiscVecT>
  void VTKOutCL::PutVector( const DiscVecT& f, const std::string& name)
/** Writes the values of a vector valued function into the VTK file
    \param name name of the function
    \param f    function*/
{
    VectorBaseCL<float> allData;
    GatherVector( f, allData);
    WriteValues(allData, name, 3);
}

template <typename DiscScalT>
  void VTKOutCL::GatherScalar(const DiscScalT& f, VectorBaseCL<float>& locData) const
/** Gathers the values of the function f on vertices and edges*/
{
    // Allocate memory
    locData.resize(numLocPoints_);

    // Get values on vertices
    Uint pos= 0;
    for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(); it!=mg_.GetTriangVertexEnd(); ++it){
        locData[pos++]= (float)f.val( *it);
    }

    // Get values on edges
    for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(); it!=mg_.GetTriangEdgeEnd(); ++it){
        locData[pos++]= (float)f.val( *it, 0.5);
    }
}

template <typename DiscVecT>
  void VTKOutCL::GatherVector(const DiscVecT& f, VectorBaseCL<float>& locData) const
/** Gather values of the function f on vertices and edges*/
{
    // Allocate memory
    locData.resize(3*numLocPoints_);

    // Get values of vertices
    Uint pos= 0;
    for (MultiGridCL::const_TriangVertexIteratorCL it= mg_.GetTriangVertexBegin(); it!=mg_.GetTriangVertexEnd(); ++it){
        const Point3DCL val= f.val( *it);
        for (int j=0; j<3; ++j)
            locData[pos++]= (float)val[j];
    }

    // Get values on edges
    for (MultiGridCL::const_TriangEdgeIteratorCL it= mg_.GetTriangEdgeBegin(); it!=mg_.GetTriangEdgeEnd(); ++it){
        const Point3DCL val= f.val( *it, 0.5);
        for (int j=0; j<3; ++j)
            locData[pos++]= (float)val[j];
    }
}

} // end of namespace DROPS

#endif
