#############################################
#   D R O P S   top level makefile          #
#############################################

# variables:

PACKAGES = geom num out misc poisson stokes navstokes tests

DROPS_ROOT = .


# rules:

all: dep $(PACKAGES:%=all_%)
	@echo "--> All executables generated successfully!"

clean: $(PACKAGES:%=clean_%)
	@echo "--> All object files and executables removed!"

distclean: $(PACKAGES:%=distclean_%)
	rm -f $(DEPFILE)
	@echo "--> Everything is cleaned up now!"

dep: topo deldepend $(PACKAGES:%=depend_%)
	@echo "--> Actual dependencies generated in $(DEPFILE)!"

check:
	cd ./tests && make check

all_%:
	cd $* && make all

clean_%:
	cd $* && make clean

distclean_%:
	cd $* && make distclean

topo:
	cd ./geom && make topo.cpp
	@echo "--> topo.cpp generated!"

deldepend:
	cp -f Dep.in $(DEPFILE)

depend_%:
	cd $* && \
        $(DEPEND) -- $(CFLAGS) -- -s"# $* dependencies:" -f- ../$*/*.cpp >> ../$(DEPFILE); \
        echo " " >> ../$(DEPFILE)

prog_%:
	cd $(@D) && make $(*F)


.PHONY: all clean distclean dep check

# include settings from the config file drops.conf:
include drops.conf
