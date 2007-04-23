/// \file
/// \brief Place for DROPS-wide parts of the docs.
///
/// The documentation main-page is in this file. Also, docs for things
/// that span multiple files:
///     - large namespaces
///     - directories
///     - bibliography
///     .


/// \namespace DROPS
/// \brief Nearly all components of DROPS are in this namespace.
///
/// The components for grid-management, discretisation, linear-algebra
/// are defined in DROPS. The source-files containing Strategy-functions
/// are not consistent in this aspect.


/// \dir ./doc
/// \brief Documentation.

/// \dir ./geom
/// \brief Multigrid-creation, -refinement; geometric boundary-description.

/// \dir ./levelset
/// \brief Components for 2-phase problems.
///
/// Programs for drop- and film- simulation reside here. Moreover, the
/// components for evolution and reparametrization of levelset-functions
/// are in the contained files.

/// \dir ./misc
/// \brief Everything, that fits in no other place.
///
/// Parameter-file handling, methods for error- and time- reporting,
/// debugging aids, etc.

/// \dir ./navstokes
/// \brief Discretisation of the Navier-Stokes-equations.

/// \dir ./num
/// \brief PDE-boundary-data-management, linear-algebra and numerical
///        quadrature.
///
/// This directory contains components for solving large, sparse systems
/// of linear equations and routines for numerical quadrature.
///
/// Also, the structures for coupling algebraic, geometric and
/// PDE-boundary data-structures is located here.

/// \dir ./out
/// Output routines

/// \dir ./poisson
/// \brief Discretisation of Poisson's-equation.

/// \dir ./stokes
/// \brief Discretisation of the Stokes-equations.


/// \mainpage The DROPS-Package for Numerical Simulations
///
/// This is the programmer's manual of the software package called
/// DROPS (cf. \ref homepage), which has recently been developed at the \ref IGPM (Institut für
/// Geometrie und Praktische Mathematik) at the Aachen University of
/// Technology. This development is still continuing and our aim is to
/// build an efficient software tool for the numerical simulation of
/// incompressible flows.
///
/// In the field of incompressible CFD already quite a few packages
/// exist. The state of the art, however, is such that for complicated
/// three-dimensional fluid dynamics (e.g. turbulent flows, multiphase
/// reacting flows, flows with free boundaries) a black-box solver is not
/// yet available. The DROPS package is developed in an interdisciplinary
/// project (SFB 540 ``Model-based Experimental Analysis of Kinetic
/// Phenomena in Fluid Multi-phase Reactive Systems'', cf. \ref SFB)
/// where complicated flow phenomena are investigated. The modeling
/// of several complex physical phenomena in this project (e.g., mass
/// transfer between liquid drops and a surrounding fluid or the fluid
/// dynamics and heat transport in a laminar falling film) requires a
/// flexible efficient and robust CFD package. Our aim is to provide
/// such a tool. From the scientific computing point of view it is of
/// interest to develop a CFD code in which several modern numerical
/// techniques which have recently been introduced in the literature are
/// implemented. Examples of such techniques are error estimation methods
/// for Navier-Stokes equations (\ref B4Eriksson, \ref B4Rannacher2), fast
/// and robust iterative solvers for discretized Navier-Stokes equations
/// with large Reynolds numbers (\ref B4Elman1, \ref B4Elman2) and level
/// set methods for two-phase flows (\ref B4Sussman, \ref B4Zhao).


/// \page bib Bibliography
/// Please find the references to the literature cited in this manual
/// below.
///
/// \section bib-links Links
/// \par SFB 540
/// \anchor SFB
/// ``Model-based Experimental Analysis of Kinetic Phenomena in Fluid
/// Multi-phase Reactive Systems'',
/// http://www.sfb540.rwth-aachen.de
///
/// \par IGPM
/// \anchor IGPM
/// Institut für Geometrie und Praktische Mathematik,
/// http://www.igpm.rwth-aachen.de
///
/// \par DROPS home page
/// \anchor homepage
/// <b>D</b>ROPS <b>ro</b>bust <b>p</b>arallel <b>s</b>olver, find some more
/// information on our home page http://www.igpm.rwth-aachen.de/DROPS
///
/// \section bib-math Mathematics
/// \par B4Eriksson \anchor B4Eriksson
/// K. Eriksson, D. Estep, P. Hansbo, and C. Johnson,
/// <em>Introduction to adaptive methods for differential equations</em>,
/// Acta Numerica 1995, pp. 105--158, Cambridge University Press, 1995.
///
/// \par B4Rannacher2 \anchor B4Rannacher2
/// R. Rannacher,
/// <em>Error control in finite element computations</em>, in: Proceedings
/// NATO-Summer School ``Error control and adaptivity in scientific computing''
/// (H. Bulgak, C. Zenger, eds.), pp. 247--278, NATO science series, series C,
/// Vol. 536. Kluwer, Dordrecht, 1999.
///
/// \par B4Elman1 \anchor B4Elman1
/// H. C. Elman and D. Silvester,
/// <em>Fast nonsymmetric iterations and preconditioning for Navier-Stokes
/// equations</em>, SIAM J. Sci. Comp. 17, pp. 33--46, 1996.
///
/// \par B4Elman2 \anchor B4Elman2
/// H. C. Elman,
/// <em>Preconditioning of the steady-state Navier-Stokes equations with low
/// viscosity</em>, SIAM J. Sci. Comp. 20, pp. 1299--1316, 1999.
///
/// \par B4Sussman \anchor B4Sussman
/// M. Sussman, P. Smereka, and S. Osher,
/// <em>A level set approach for computing solutions to incompressible two-phase
/// flow</em>, J. Comp. Phys. 114, pp. 146--159, 1994.
///
/// \par B4Zhao \anchor B4Zhao
/// H.-K. Zhao, B. Merriman, S. Osher, and L. Wang,
/// <em>Capturing the behaviour of bubbles and drops using the variational level
/// set approach</em>, J. Comp. Phys. 143, pp. 495--518, 1998.

