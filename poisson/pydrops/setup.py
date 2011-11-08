from distutils.core import setup, Extension

numpy_include = ""
import numpy; numpy_include = numpy.get_include()


drops_module = Extension('drops', sources=['py_interpolation.cpp',
                                           'py_utils.cpp',
                                           'pydrops.cpp'],
                         include_dirs=[numpy_include,
                                       '../../'],
#                                       '../../../poisson',
#                                       '../../../geom',
#                                      '../../../num',
#                                       '../../../misc'],
                         extra_objects=['../../geom/topo.o',
                                        '../../geom/multigrid.o',
                                        '../../geom/boundary.o',
                                        '../../geom/builder.o',
                                        '../../num/discretize.o',
                                        '../../num/unknowns.o',
                                        '../../misc/problem.o',
                                        '../../misc/utils.o',
                                        '../../poisson.o',
                                        'drops_utils.o',
                                        'py_coeff_dp_stat.o',
                                        'py_source.o',
                                        'py_kdelta_psi.o'
                                        ],
                         depends=['py_source.hpp',
                                  'py_source.cpp',
                                  'py_param.hpp',
                                  'py_utils.hpp',
                                  'py_interpolation.hpp',
                                  'py_L2_alt_scalar_prod.cpp',
                                  'functors.hpp',
                                  'py_functors.hpp',
                                  'py_coeff_gradient.hpp',
                                  'py_coeff_gradient.cpp',
                                  'py_coeff_dp_stat.hpp',
                                  'drops_utils.hpp',
                                  '../../poisson.h',
                                  '../../poisson.cpp',
                                  '../../instatpoisson.tpp',
                                  '../../integrTime.h',
                                  '../../../num/solver.h'],
                         extra_compile_args=['-W','-Wall','-pedantic', '-O0', '-g', '-fopenmp'])

setup(name='drops', version='1.0', ext_modules=[drops_module])
