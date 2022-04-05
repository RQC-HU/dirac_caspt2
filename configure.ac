AC_PREREQ([2.69])
AC_INIT([relcaspt2_parallel],[0.4.5],[kohei-noda@hiroshima-u.ac.jp])
AC_CONFIG_SRCDIR([src])

#デフォルトのコンパイルオプション
DEFAULT_FCFLAGS="-cpp "

#FCFLAGS="-xHost -mkl -cpp -DDIRAC21 -pad -Ofast -mp1 -integer_size 64 -unroll -qopenmp -warn all -warn nounused -nogen-interfaces"
FCFLAGS_PRESERVE=$FCFLAGS
#使用するコンパイラの決定(./configure FC=ifort などとして使用者がコンパイラを指定することもできる)
AC_PROG_FC
AC_FC_MODULE_EXTENSION
AC_FC_MODULE_FLAG
AC_FC_MODULE_OUTPUT_FLAG

#コンパイラの実行ファイルの名前の取得
FCNAME=`echo $FC | rev | cut -d "/" -f 1 | rev`

#Intel Fortranのデフォルトコンパイルオプション
INTEL_FORTRAN_FLAGS='-mkl -pad -mp1 -integer_size 64 -unroll -qopenmp -warn all -warn nounused -nogen-interfaces'
OPT='-Ofast '
#Mklのテスト

AC_CHECK_HEADERS([mkl.h], [echo 'Found mkl.h!' ], [AC_MSG_ERROR([mkl.h not found or not working])], [])
#gnu Fortran用のmklリンク(少なくともrelqc01では使用可能)
LDFLAGS='-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm'
#LDFLAGS=`-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_ilp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl`
AS_IF(
	[test $FCNAME == gfortran],
	[DEFAULT_FCFLAGS+='-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm -fopenmp -fdefault-integer-8 -m64'],
	[test $FCNAME == mpiifort],
	[DEFAULT_FCFLAGS+=$OPT"-DHAVE_MPI "$INTEL_FORTRAN_FLAGS],
	[test $FCNAME == mpifort],
    [DEFAULT_FCFLAGS+=$OPT'-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm -DHAVE_MPI -fopenmp']
	[test $FCNAME == ifort],
	[DEFAULT_FCFLAGS+=$OPT$INTEL_FORTRAN_FLAGS]
	[echo 'Your fortran compiler is not gfortran or ifort']
)
AC_SUBST([DEFAULT_FCFLAGS])
AM_INIT_AUTOMAKE([foreign subdir-objects])

AC_CONFIG_FILES([Makefile
                 src/Makefile])

#Option for DIRAC version >= 21
AC_ARG_ENABLE([mpi],
	[AS_HELP_STRING([--enable-mpi=yes|no],
	[If you use mpi, Enable this feature (default is no)])],
	[],
	[enable_mpi="no"]
)

#Option for DIRAC version >= 21
AC_ARG_ENABLE([dirac21],
	[AS_HELP_STRING([--enable-dirac21=yes|no],
	[If you use DIRAC version >= 21, Enable this feature (default is no)])],
	[],
	[enable_dirac21="no"]
)
AS_IF([test "x$enable_dirac21" = xyes],
	[DIRAC21=-DDIRAC21],
	[DIRAC21=]
)
if test "${FCFLAGS}" == '-g'; then
FCFLAGS=''
fi

echo $FCFLAGS

AC_MSG_CHECKING(enable_dirac21)
AC_SUBST([DIRAC21])
AC_OUTPUT
