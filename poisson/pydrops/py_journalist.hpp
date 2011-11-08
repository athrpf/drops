
#ifndef __PY_JOURNALIST_HPP__
#define __PY_JOURNALIST_HPP__

#include <fstream>
#include <iostream>
#include <string>

enum PrintLevel {
  JNLST_PRINTLEVEL_1=1,    // print only the most important info
  JNLST_PRINTLEVEL_ALL=-1, // print everything
  JNLST_PRINTLEVEL_NONE=0  // print nothing
};

class Journalist
{
private:
  PrintLevel level_;
  Journalist();
public:
  Journalist(PrintLevel level=JNLST_PRINTLEVEL_1) : level_(level) {}


  friend Journalist& operator<<(Journalist& jnlst, std::ostream& out){
    if (jnlst.level_!=0) {
      std::cout << out;
    }
    return jnlst;
  }
  friend Journalist& operator<<(Journalist& jnlst, const std::string& out){
    if (jnlst.level_!=0) {
      std::cout << out;
    }
    return jnlst;
  }
  friend Journalist& operator<<(Journalist& jnlst, const double& out){
    if (jnlst.level_!=0) {
      std::cout << out;
    }
    return jnlst;
  }
  friend Journalist& operator<<(Journalist& jnlst, const int& out){
    if (jnlst.level_!=0) {
      std::cout << out;
    }
    return jnlst;
  }

};


#endif
