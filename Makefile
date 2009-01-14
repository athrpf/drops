#############################################
#   D R O P S   top-level makefile          #
#############################################

# variables:

PARPACKAGES = #parallel partests
PACKAGES = geom num out misc poisson stokes navstokes tests levelset surfactant $(PARPACKAGES)

DROPS_ROOT = .


# rules:

default: dep all

#all: $(PARPACKAGES:%=all_%)
all: $(PACKAGES:%=all_%)
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


.PHONY: all clean distclean distclean_dox default dep deldepend doc stat topo check

# include settings from the config file drops.conf:
include drops.conf

DDD:
	cd $(DDD_HOME) && ./install && gmake clean && \
	gmake -j ARCH_CFLAGS="$(OPTFLAGS)" ARCH_TYPE="__PC__" ARCH_CC="$(ARCH_CC)" ARCH_LINK="$(ARCH_CC)" && \
	cd $(DROPS_ROOT)

ParMetis:
	cd $(PARMETIS_HOME) && gmake clean && \
	gmake COPTIONS="$(OPTFLAGS)" CC="$(ARCH_CC)" LD="$(ARCH_CC)" && cd $(DROPS_ROOT)
