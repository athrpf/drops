#
# cxxopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  The version tested against
#               is specified by platform.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Copyright 1984-2000 The MathWorks, Inc.
# $Revision$  $Date$
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex"
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The cmex script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        alpha)
#----------------------------------------------------------------------------
#           cxx -V
#           Compaq C++ V6.2-024 for Digital UNIX V4.0F  (Rev. 1229)
            CC='cxx'
            CFLAGS='-shared -ieee -pthread'
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,-expect_unresolved,'*',-hidden,-exported_symbol,$ENTRYPOINT,-exported_symbol,mexVersion,-exported_symbol,'__*'"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        hpux)
#----------------------------------------------------------------------------
#           what `which aCC`
#           HP aC++ B3910B A.03.13
#           HP aC++ B3910B X.03.11.10 Language Support Library
            CC='aCC'
            CFLAGS='+Z +DA2.0 -D_POSIX_C_SOURCE=199506L -D_HPUX_SOURCE'
            CLIBS="$MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-b -Wl,+e,$ENTRYPOINT,+e,mexVersion,+e,_shlInit"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        hp700)
#----------------------------------------------------------------------------
#           what `which aCC`
#           HP aC++ B3910B A.01.21
#           HP aC++ B3910B A.01.19.02 Language Support Library
            CC='aCC'
#           Remove +DAportable from CFLAGS if you wish to optimize
#           for target machine
            CFLAGS='+Z -D_HPUX_SOURCE +DAportable'
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-b -Wl,+e,$ENTRYPOINT,+e,mexVersion,+e,_shlInit,+e,errno"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        ibm_rs)
#----------------------------------------------------------------------------
#           lslpp -l ibmcxx.cmp
#           3.6.6.0  COMMITTED  IBM C and C++ Compilers
            CC='/usr/vacpp/bin/xlC'
            CFLAGS='-D_THREAD_SAFE -D_ALL_SOURCE -qrtti=all'
            CLIBS="$MLIBS -lmatlb -lmat -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD='/usr/vacpp/bin/makeC++SharedLib'
            LDFLAGS='-p 0'
            LDOPTIMFLAGS=''
            LDDEBUGFLAGS=''
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,--rpath-link,$TMW_ROOT/extern/lib/$Arch,--rpath-link,$TMW_ROOT/bin/$Arch"
#           g++ -v
#           gcc version 2.95.2 19991024 (release)
            CC='g++'
#           Add -fhandle-exceptions to CFLAGS to support exception handling
            CFLAGS='-fPIC -ansi -D_GNU_SOURCE -pthread'
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sgi)
#----------------------------------------------------------------------------
#           CC -version
#           MIPSpro Compilers: Version 7.3.1.2m
            CC='CC'
#           Add -exceptions to CFLAGS to support exception handling
            CFLAGS='-n32 -D_POSIX_C_SOURCE=199506L -D__EXTENSIONS__ -D_XOPEN_SOURCE -D_XOPEN_SOURCE_EXTENDED'
            CLIBS="-dont_warn_unused $MLIBS -lm -lC"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-n32 -shared -exported_symbol $ENTRYPOINT -exported_symbol mexVersion"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol2)
#----------------------------------------------------------------------------
#           CC -V
#           WorkShop Compilers 5.0 98/12/15 C++ 5.0
            CC='CC -compat=5'
            CCV=`CC -V 2>&1`
            version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
            if [ "$version" = "4" ]; then
                    echo "SC5.0 or later C++ compiler is required"
            fi
            CFLAGS='-KPIC -dalign -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt'
            CLIBS="$MLIBS -lm -lCstd -lCrun"

            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDFLAGS="-G -mt -M$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
