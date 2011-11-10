#############################################
#   D R O P S   top-level makefile          #
#############################################
DROPS_ROOT = .

# include settings from the config file drops.conf:
include drops.conf

# variables:

PARPACKAGES = parallel partests poisson levelset
SERPACKAGES = geom num out misc poisson stokes navstokes tests levelset surfactant transport
PACKAGES = $(SERPACKAGES) $(PARPACKAGES)
BUILDPACKAGES = $(if $(PAR_BUILD),$(PARPACKAGES),$(SERPACKAGES))

# rules:
pydrops: pydrops_deps
	cd poisson/pydrops; python setup.py build

pydrops_deps:
	cd misc; $(MAKE) problem.o utils.o
	cd geom; $(MAKE) topo.o multigrid.o boundary.o builder.o
	cd num; $(MAKE) discretize.o unknowns.o
	cd out; $(MAKE) output.o
	cd poisson; $(MAKE) poisson.o
	cd poisson/pydrops; $(MAKE) drops_utils.o py_coeff_dp_stat.o py_source.o py_kdelta_psi.o

default: dep all

all: $(BUILDPACKAGES:%=all_%)
	@echo "--> All executables generated successfully!"

strip: $(PACKAGES:%=strip_%)

clean: $(PACKAGES:%=clean_%)
	@echo "--> All object files and executables removed!"

distclean: $(PACKAGES:%=distclean_%) distclean_dox
	rm -f $(DEPFILE)
	@echo "--> Everything is cleaned up now!"

dep: deldepend $(PACKAGES:%=depend_%)
	@echo "--> Actual dependencies generated in $(DEPFILE)!"

check:
	cd ./tests && $(MAKE) check

doc:
	doxygen dox.cfg

stat:
	-ls $(PACKAGES:%=%/*.[ct]pp) $(PACKAGES:%=%/*.h) > $(DROPS_ROOT)/srcFiles.tmp
	cat $(DROPS_ROOT)/srcFiles.tmp | xargs wc --lines
	@rm -f $(DROPS_ROOT)/srcFiles.tmp
	@echo "--> lines of code (LOC) statistics done"

all_%:
	cd $* && $(MAKE) all

strip_%:
	cd $* && $(MAKE) strip

clean_%:
	cd $* && $(MAKE) clean

distclean_%:
	cd $* && $(MAKE) distclean

distclean_dox:
	cd ./doc && rm -rf dox

clean_DDD:
	cd $(DDD_HOME) && gmake clean && cd $(DROPS_ROOT)

clean_ParMetis:
	cd $(PARMETIS_HOME) && gmake clean && cd $(DROPS_ROOT)

clean_HYPRE:
	cd $(HYPRE_HOME) && gmake clean && cd $(DROPS_ROOT)

topo:
	cd ./geom && $(MAKE) topo.cpp
	@echo "--> topo.cpp generated!"

deldepend:
	cp -f Dep.in $(DEPFILE)

depend_%:
	cd $* && \
        $(DEPEND) -- $(CXXFLAGS) -- -s"# $* dependencies:" -f- ../$*/*.cpp >> ../$(DEPFILE) 2>/dev/null; \
        echo " " >> ../$(DEPFILE)

prog_%:
	cd $(@D) && $(MAKE) $(*F)

libs: DDD ParMetis

clean_libs: clean_DDD clean_ParMetis

DDD:
	cd $(DDD_HOME) && ./install && gmake clean && \
	gmake -j ARCH_CFLAGS="$(OPTFLAGS)" ARCH_TYPE="__PC__" ARCH_CC="$(ARCH_CC)" ARCH_LINK="$(ARCH_CC)"

ParMetis:
	cd $(PARMETIS_HOME) && gmake clean && \
	gmake COPTIONS="$(OPTFLAGS)" CC="$(ARCH_CC)" LD="$(ARCH_CC)"

HYPRE:
	cd $(HYPRE_HOME) && gmake install

.PHONY: all clean distclean distclean_dox default dep deldepend doc stat topo check

