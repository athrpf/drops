

#include "py_utils.hpp"

bool check_dimensions_3(PyArrayObject* obj, int N1, int N2, int N3, const std::string name)
{
  if (obj->nd!=3){
    std::stringstream ss;
    ss << name << " has " << obj->nd << " dimensions, but should have 3!\n";
    PyErr_SetString(PyExc_ValueError, ss.str().c_str());
    return false;
  }
  npy_intp* od = obj->dimensions;
  if (!(od[0]==N1 && od[1]==N2 && od[2]==N3)) {
    std::stringstream ss;
    ss << name << " has dimensions (" << od[0] << "," << od[1] << "," << od[2]
       << ") but should be (" << N1 << "," << N2 << "," << N3 << ")\n";
    PyErr_SetString(PyExc_ValueError, ss.str().c_str());
    return false;
  }
  return true;
}

bool check_dimensions_4(PyArrayObject* obj, int N1, int N2, int N3, int N4, const std::string name)
{
  if (obj->nd!=4){
    std::stringstream ss;
    ss << name << " has " << obj->nd << " dimensions, but should have 4!\n";
    PyErr_SetString(PyExc_ValueError, ss.str().c_str());
    return false;
  }
  npy_intp* od = obj->dimensions;
  if (!(od[0]==N1 && od[1]==N2 && od[2]==N3 && od[3]==N4)) {
    std::stringstream ss;
    ss << name << " has dimensions (" << od[0] << "," << od[1] << "," << od[2] << ","
       << od[3]<<") but should be (" << N1 << "," << N2 << "," << N3 << "," << N4 << ")\n";
    PyErr_SetString(PyExc_ValueError, ss.str().c_str());
    return false;
  }
  return true;
}
