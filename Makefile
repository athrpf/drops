#############################################
#   D R O P S   top level makefile          #
#############################################

# variables:

PACKAGES = geom num out misc poisson stokes navstokes tests


# rules:

all: dep $(PACKAGES:%=all_%)
	@echo "--> All executables generated succesfully!"
        
clean: $(PACKAGES:%=clean_%)
	@echo "--> All object files and executables removed!"
        
distclean: $(PACKAGES:%=distclean_%)
	rm -f $(DEPFILE)
	@echo "--> Everything is cleaned up now!"
        
dep: deldepend $(PACKAGES:%=depend_%)
	@echo "--> Actual dependencies generated in $(DEPFILE)!"
        
all_%:
	cd $*; make all

clean_%:
	cd $*; make clean
        
distclean_%:
	cd $*; make distclean
        
deldepend:
	cp -f Dep.in $(DEPFILE)

depend_%:
	cd $*; \
        $(DEPEND) -- $(CFLAGS) -- -s"# $* dependencies:" -f- ../$*/*.cpp >> ../$(DEPFILE); \
        echo " " >> ../$(DEPFILE)
        
prog_%:
	cd $(@D); make $(*F)

        
.PHONY: all clean distclean dep

# include settings from the config file drops.conf:
DROPS_ROOT = .
include drops.conf
