/***************************************************************************
*  File:    compensight.cpp                                                *
*  Content: Comparing ensight-files                                        *
*  Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
*           Oliver Fortmeier, RZ RWTH Aachen                               *
*  begin:           20.10.2006                                             *
*  last modified:   30.10.2006                                             *
*                   27.11.2006 (relative diff)                             *
***************************************************************************/
/// \author Oliver Fortmeier
/// \file compensight.cpp
/// \brief Comparing ensight-files

#include "misc/utils.h"
#include "misc/container.h"
#include "misc/params.h"

#include <map>
#include <list>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

namespace DROPS
{

/// \brief Lexicographic order
/** This specialized function uses lexicographical odering for comparing two SVectorCL<3>.*/
template<>
  inline bool operator< (const SVectorCL<3>& a, const SVectorCL<3>& b)
{
    for (Uint i=0; i<3; ++i)
        if (a[i] < b[i])
            return true;
        else if ( a[i] > b[i])
            return false;
    return false;
}

/// \brief Read a point from an input stream
template<Uint size>
  inline std::istream& operator>> (std::istream& is, SVectorCL<size>& a)
{
    for (Uint i=0; i<size; ++i)
        is >> (a[i]);
    return is;
}

/// \brief Class that reads ensight files and compares them
template<typename T>
class ValueCL
/// This class can read ensight geometry- and value files. Also
/// one can obtain values on a given node. The template parameter
/// spezifies if the value on a node is scalar or a vector.
///
/// This class also remember, in which ordering the values are
/// read. Hence this programm can also be used for transforming
/// an ensight data file to another geometry file.
/// \todo (of): Using param-Files for specification of the input files, scalar or vectoriell mode, binary mode ...
{
  public:
    /// \brief type of the value on a node
    typedef T EntryT;
    /// \brief The values are stored within a map which key is the node and value is the value on the node
    typedef std::map<Point3DCL,EntryT>       ValueCT;
    /// \brief A single element of the map
    typedef std::pair<Point3DCL,EntryT>      ValueT;
    /// \brief constant iterator over map
    typedef typename ValueCT::const_iterator ValueIT;
    /// \brief Storing the ordering of the input, is done in a list
    typedef std::list<ValueIT>               ListCT;
    /// \brief constant iterator over the orderer list
    typedef typename ListCT::const_iterator  ListIT;

  private:
    mutable ValueCT vals_;
    ListCT orderList_;
    std::ifstream geomFile_;
    std::ifstream valFile_;

  public:
    ValueCL(const std::string&, const std::string&);            ///< Constructor with given files

    template<typename FileT>
    bool CheckFile(const FileT&) const;                         ///< Check if a file is correct opend
    size_t Size() const { return vals_.size(); }                ///< Get number of elements

    void Read();                                                ///< Read geometry and values

    const EntryT& operator[] (const Point3DCL& p) const         ///< Get value on a point
        { return vals_[p]; }

    ValueIT begin() const { return vals_.begin(); }             ///< first element
    ValueIT end()   const { return vals_.end(); }               ///< behind last element

    ListIT listBegin() const { return orderList_.begin(); }     ///< first element of ordered list
    ListIT listEnd() const { return orderList_.end(); }         ///< end of the ordered list
};

template<typename T>
  ValueCL<T>::ValueCL(const std::string& geomFileName, const std::string& valFileName)
    : geomFile_(geomFileName.c_str()), valFile_(valFileName.c_str())
/// \param geomFileName name of the file that contains geometry information
/// \param valFileName  name of the file that contains value information
{
    if (!CheckFile(geomFile_)){
        std::cout << "Cannot open file "<<geomFileName<<std::endl;
        exit(-1);
    }

    if (!CheckFile(valFile_)){
        std::cout << "Cannot open file "<<valFileName<<std::endl;
        exit(-1);
    }
}

/// \brief Check if the file could be read
template<typename T> template<typename FileT>
  bool ValueCL<T>::CheckFile(const FileT& is) const
/// Check if the file exists and could be opened
/// \param is file that should be checked
{
    return !(!is);
}

template<typename T>
  void ValueCL<T>::Read()
/// Just read all values that are in the given file
/// \pre Right now the geometry file must not contain parts
{
    const int describerLinesGeom = 5;
    const int describerLinesValue= 1;
    char buf[256];

    for (int i=0; i<describerLinesGeom; ++i)     // ignore describer lines of geometry file
        geomFile_.getline( buf, 256);

    for (int i=0; i<describerLinesValue; ++i)    // ignore describer lines of value file
        valFile_.getline( buf, 256);

    int numNodes=-1;
    int dummy;
    geomFile_ >> numNodes;

    Point3DCL node;
    T         val;
    for (int i=0; i<numNodes; ++i){
        geomFile_ >> dummy;
        geomFile_ >> node;
        valFile_  >> val;
        vals_[node] = val;
        orderList_.push_back(vals_.find(node));
    }
}

/// \brief Get scalar value on a node
template<typename T>
  double getScalarVal(const T& a);

/// \brief Get scalar value on a node with scalar values
template<>
  inline double getScalarVal(const double& a)
/// just returns the value
{
    return a;
}

/// \brief Get scalar value on a node with vectoriell values
template<>
  inline double getScalarVal(const Point3DCL& a)
/// returns squared euklidian norm of the 3D vector
{
    return a.norm_sq();
}

/// \brief distance between two values
template<typename T>
  double my_distance(const T& a, const T& b)
{
    return std::fabs(getScalarVal(a-b));
}

/// \brief maximal value
template<typename T>
  double my_max(const T& a, const T& b)
{
    return std::max(getScalarVal(a), getScalarVal(b));
}

/// \brief minimal value
template<typename T>
  double my_min(const T& a, const T& b)
{
    return std::min(getScalarVal(a), getScalarVal(b));
}

/// \brief relative difference of two values
template<typename T>
  double reldiff(const T& a, const T& b)
/// relative difference is difference divided by bigger value
{
    const double val_a=getScalarVal(a);
    const double val_b=getScalarVal(b);
    return std::abs(val_a-val_b)/std::max(std::abs(val_a),std::abs(val_b));
}

/// \brief Check differences between two ValueCL's
template<typename T>
  void Compare (const ValueCL<T>& a, const ValueCL<T>& b)
/** This function compares two ValueCL. Therefore the following aspects are considered:
    <ul>
        <li>biggest absolute difference on a node</li>
        <li>biggest relative difference on a node</li>
        <li>euklidian norm between whole vector</li>
        <li>weighted euklidian norm between whole vector: sqrt(sum_i=0^(n-1) ||a(i)-b(i)||^2/n)</li>
        <li>max value</li>
        <li>min value</li>
        <li>number of nodes</li>
    </ul>
   \pre both ValueCL must have same length
*/
{
    if (a.Size()!=b.Size()){
        std::cout << "Not the same length!\n";
        return;
    }

    double diff, rel_diff, max_diff=-1., rel_max_diff=-1.,      // differences
           norm2=0, normw2=0,                                   // norms
           max_val=0, min_val=1e99;                             // max and min val
    T val1=T(), val2=T();                                       // values of max diff
    T rel_val1=T(), rel_val2=T();                               // values of max_rel diff
    Point3DCL max_node, rel_max_node;                           // nodes of max diff and max_rel diff

    for (typename ValueCL<T>::ValueIT it(a.begin()), end(a.end()); it!=end; ++it)
    {   // compare all values

        // Get values
        T first=it->second,
          second=b[it->first];

        // differences
        diff = my_distance(first, second);
        rel_diff=reldiff(first, second);

        if (diff>max_diff)
        {   // new biggest max diff
            max_diff=diff;
            max_node=it->first;
            val1= it->second;
            val2= b[it->first];
        }

        if (rel_diff>rel_max_diff)
        {   // new biggest max_rel diff
            rel_max_diff=rel_diff;
            rel_max_node=it->first;
            rel_val1= it->second;
            rel_val2= b[it->first];
        }

        norm2 += diff*diff;                                 // euklidian norm
        max_val= std::max(max_val, my_max(first, second));
        min_val= std::min(min_val, my_min(first, second));
    }

    normw2= std::sqrt(norm2/a.Size());
    norm2 = std::sqrt(norm2);

    // output
    std::cout << " Summary:"
              << "\n max_diff          = "<<max_diff<<" on node "<<max_node<<": "<<val1<<" instead of "<<val2
              << "\n relative max_diff = "<<rel_max_diff<<" on node "<<rel_max_node<<": "<<rel_val1<<" instead of "<<rel_val2
              << "\n 2-norm            = "<<norm2
              << "\n w2-norm           = "<<normw2
              << "\n max_value         = "<<max_val
              << "\n min_value         = "<<min_val<<"\n"
              << "\n number of values  = "<<a.Size()
              << std::endl;
}

/// \brief Write a Point or a double in ensight format onto an output stream
template<typename T>
  void write(const T&, std::ostream&, Uint& count);

template<>
  inline void write(const double& x, std::ostream& os, Uint& count)
/** Writes a double onto an outputstream and increment counter by one
    \param x     variable for output
    \param os    output stream
    \param count counter
*/
{
    os << std::setw(12) << x;
    ++count;
}


template<>
  inline void write(const Point3DCL& p, std::ostream& os, Uint& count)
/** Writes vectoriel data onto an outputstream and increment counter by three
    \param p     variable for output
    \param os    output stream
    \param count counter
*/
{
    for (int i=0; i<3; ++i)
        os << std::setw(12) << p[i];
    count+=3;
}


/// \brief Transform one ensight data file according to the gemetry of another ensight file
template<typename T>
  void Transform (const ValueCL<T>& refVal, const ValueCL<T>& toTransVal, const std::string filename)
/**
    This function takes the values out of \a toTransVal and writes them in the ordering, so that is a
    valid ensight data-file for the geometry of refVal.
    \param refVal     ValueCL where the geometry is stored, which ensight can visualise
    \param toTransVal Data file, that should be transformed to a valid ensight data file
    \param filename   Name of the file where to store the ensight data file
*/
{
    if (refVal.Size()!=toTransVal.Size()){
        std::cout << "Transform: Incompatible sizes"<<std::endl;
        std::exit(-1);
    }

    std::ofstream out(filename.c_str());
    refVal.CheckFile(out);
    out.flags(std::ios_base::scientific);
    out.precision(5);
    out.width(12);

    out << "DROPS data file:" << std::endl;

    Uint cnt=0;
    for (typename ValueCL<T>::ListIT it(refVal.listBegin()); it!=refVal.listEnd(); ++it)
    {
        write(toTransVal[(**it).first], out, cnt);
        if ( cnt==6 )
        { // Ensight expects six real numbers per line
            cnt= 0;
            out << '\n';
        }
    }
    out << std::endl;
    out.close();
}

} // end of namespace DROPS

using DROPS::ValueCL;
using DROPS::Compare;
using DROPS::Transform;

int main(int argc, char* argv[])
{
    if (argc!=6 && argc!=8){
        std::cout << "Usage (compare): "<<argv[0]<<" <0 scalar/ 1 vector> <par_geom_file> <par_val_file> "
                  << "<ser_geom_file> <ser_val_file> "
                  <<std::endl<<"or"<<std::endl
                  << "Usage(transform): "<<argv[0]<<" -t <0 scalar/ 1 vector> <par_geom_file> <par_val_file> "
                  << "<ser_geom_file> <ser_val_file> <out_val_file>"
                  <<std::endl;
        exit(-1);
    }

    std::string mode(argv[1]);
    int start_idx=(mode=="-t" ? 2 : 1);

    int type=std::atoi(argv[start_idx]);
    std::string a(argv[start_idx+1]), b(argv[start_idx+2]),
                c(argv[start_idx+3]), d(argv[start_idx+4]);

    std::cout << "Used geometry:"<<a<<" and "<<c<<std::endl;
    std::cout << "Used data:    "<<b<<" and "<<d<<std::endl;

    if (type)    // vectoriel data
    {
        std::cout << "Using vectoriel datas\n";
        ValueCL<DROPS::Point3DCL> ParVals(a,b);
        ValueCL<DROPS::Point3DCL> SerVals(c,d);
        std::cout << "Reading parallel values ... \n";
        ParVals.Read();
        std::cout << "Reading serial values ... \n";
        SerVals.Read();
        if (mode!="-t")
        {
            std::cout << "Comparing values ... \n";
            Compare(ParVals,SerVals);
        }
        else
        {
            std::cout << "Transforming values ... \n";
            const std::string out_name(argv[7]);
            Transform<DROPS::Point3DCL>(ParVals, SerVals, out_name);
        }
    }
    else            // scalar data
    {
        std::cout << "Using scalar datas\n";
        ValueCL<double> ParVals(a,b);
        ValueCL<double> SerVals(c,d);
        std::cout << "Reading parallel values ... \n";
        ParVals.Read();
        std::cout << "Reading serial values ... \n";
        SerVals.Read();
        if (mode!="-t")
        {
            std::cout << "Comparing values ... \n";
            Compare(ParVals,SerVals);
        }
        else
        {
            std::cout << "Transforming values ... \n";
            const std::string out_name(argv[7]);
            Transform<double>(ParVals, SerVals, out_name);
        }

    }

    return 0;
}
