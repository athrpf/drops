#ifndef __PYPDEFUNCTION_HPP__
#define __PYPDEFUNCTION_HPP__

class PyPdeFunction : public PdeFunction {
private:
  int nx, ny, nz, nt;
  boost::python::numeric::array* data;
public:

  PyPdeFunction(boost::python::numeric::array* data_) : data(data_)
  {
    using namespace boost::python;
    tuple t = extract<tuple>(data->attr("shape"));
    assert(len(t)==4);
    nx = boost::python::extract<int>(t[0]);
    ny = boost::python::extract<int>(t[1]);
    nz = boost::python::extract<int>(t[2]);
    nt = boost::python::extract<int>(t[3]);
  }


  virtual ~PyPdeFunction(){}

  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const {
    bool retval = false;
    if (Nx==nx && Ny==ny && Nz==nz && Nt==nt) {
      retval = true;
    }
    Nx = nx; Ny = ny; Nz = nz; Nt = nt;
    return retval;
  }

  virtual double operator()(int ix, int iy, int iz, int it) const {
    assert (ix>=0 && ix<nx);
    assert (iy>=0 && iy<ny);
    assert (iz>=0 && iz<nz);
    assert (it>=0 && it<nt);
    //std::cout << "reading domain function " << ix << ","  << iy << ","  << iz << ","  << it << std::endl;
    using namespace boost::python;
    tuple t = make_tuple<int,int,int,int>(ix,iy,iz,it);
    return extract<double>((*data)[t]);
  }

  double at(int ix, int iy, int iz, int it) const {
    assert (ix>=0 && ix<nx);
    assert (iy>=0 && iy<ny);
    assert (iz>=0 && iz<nz);
    assert (it>=0 && it<nt);
    using namespace boost::python;
    tuple t = make_tuple<int,int,int,int>(ix,iy,iz,it);
    return extract<double>((*data)[t]);
  }
};

class PyPdeBoundaryFunction : public PdeFunction {
private:
  int dead_dim;
  int nx, ny, nz, nt;
  boost::python::numeric::array* data;
public:
  PyPdeBoundaryFunction(boost::python::numeric::array* data_, int dead_dim_) : dead_dim(dead_dim_), data(data_)
  {
    using namespace boost::python;
    //std::cout << "In PyPdeboundaryfunction" << std::endl;
    tuple t = extract<tuple>(data->attr("shape"));
    if (!(len(t)==3)) {
      std::cout << "Length of Boundary function array is not 3 but " << len(t) << std::endl;
      assert(false);
    }
    if (dead_dim==0) {
      nx = 1;
      ny = boost::python::extract<int>(t[0]);
      nz = boost::python::extract<int>(t[1]);
      nt = boost::python::extract<int>(t[2]);
    } else if (dead_dim==1) {
      nx = boost::python::extract<int>(t[0]);
      ny = 1;
      nz = boost::python::extract<int>(t[1]);
      nt = boost::python::extract<int>(t[2]);
    } else if (dead_dim==2) {
      nx = boost::python::extract<int>(t[0]);
      ny = boost::python::extract<int>(t[1]);
      nz = 1;
      nt = boost::python::extract<int>(t[2]);
    } else if (dead_dim==3) {
      nx = boost::python::extract<int>(t[0]);
      ny = boost::python::extract<int>(t[1]);
      nz = boost::python::extract<int>(t[2]);
      nt = 1;
    }
  }

  virtual ~PyPdeBoundaryFunction(){}

  virtual bool get_dimensions(int& Nx, int& Ny, int& Nz, int& Nt) const {
    bool retval = false;
    if (Nx==nx && Ny==ny && Nz==nz && Nt==nt) {
      retval = true;
    }
    Nx = nx; Ny = ny; Nz = nz; Nt = nt;
    return retval;
  }

  virtual double operator()(int ix, int iy, int iz, int it) const {
    return this->at(ix,iy,iz,it);
  }

  double at(int ix, int iy, int iz, int it) const {
    //std::cout << "reading boundary function " << ix << ","  << iy << ","  << iz << ","  << it << std::endl;
    //std::cout << "size is " << nx << ","  << ny << ","  << nz << ","  << nt << std::endl;
    using namespace boost::python;
    tuple t;
    if (dead_dim==0) {
      assert(iy>=0 && iy<ny);
      assert(iz>=0 && iz<nz);
      assert(it>=0 && it<nt);
      t = make_tuple<int,int,int>(iy,iz,it);
    } else if(dead_dim==1) {
      assert (ix>=0 && ix<nx);
      assert (iz>=0 && iz<nz);
      assert (it>=0 && it<nt);
      t = make_tuple<int,int,int>(ix,iz,it);
    } else if(dead_dim==2) {
      assert (ix>=0 && ix<nx);
      assert (iy>=0 && iy<ny);
      assert (it>=0 && it<nt);
      t = make_tuple<int,int,int>(ix,iy,it);
    } else if(dead_dim==3) {
      //std::cout << "evaluating initial value\n";
      assert (ix>=0 && ix<nx);
      assert (iy>=0 && iy<ny);
      assert (iz>=0 && iz<nz);
      t = make_tuple<int,int,int>(ix,iy,iz);
    }
    return extract<double>((*data)[t]);
  }
};

#endif
