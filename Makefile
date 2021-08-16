# include  /home/minori/PROGRAMS/utchem_rq37_new/utchem/config/makeconfig

#OBJS4 = ../netlibfiles/linpack/linpack.a
#R4DCASCI = four_caspt2_module.o nbitsa.o readvec.o read1mo.o \
	readorb_enesym.o \
	one_e_exct.o dim.o diag.o cutoff.o \
	readint2_casci.o prtoutfock.o fockhf1.o fockcasci.o \
	tramo.o e0after_tra.o trac.o sdet.o \
	calce0.o fockdiag.o \
	takekr.o \
	intmo.o \
	timing.o mem.o \
	uramda_s_half.o nrintread.o \
	checkdgc.o e0test_v2.o casci.o casdet.o casmat.o r4dcasci.o \
	create_binmdcint.o

#R4DCASCI_TY = four_caspt2_module.o nbitsa.o readvec.o read1mo_ty.o \
	readorb_enesym_ty.o \
	one_e_exct.o dim.o diag.o cutoff.o \
	readint2_casci_co.o prtoutfock.o fockhf1_ty.o fockcasci_ty.o \
	tramo_ty.o e0after_tra_ty.o trac.o sdet.o \
	calce0.o fockdiag_ty.o \
	takekr.o \
	intmo_ty.o \
	pgsym_ty.o \
	timing.o mem.o \
	uramda_s_half.o nrintread.o \
	checkdgc.o e0test_v2.o casci_ty.o casdet_ty.o casmat.o r4dcasci_ty.o

R4DCASCI_CO = four_caspt2_module.o nbitsa.o readvec.o read1mo_co.o \
	readorb_enesym_co.f \
	one_e_exct.o dim.o diag.o cutoff.o \
	readint2_casci_co.o prtoutfock.o fockhf1_ty.o fockcasci_ty.o \
	tramo_ty.o e0after_tra_ty.o trac.o sdet.o \
	calce0.o fockdiag_ty.o \
	takekr.o \
	intmo_ty.o \
	timing.o mem.o \
	uramda_s_half.o nrintread.o \
	checkdgc.o e0test_v2.o casci_ty.o casdet_ty.o casmat.o r4dcasci_co.o \
	create_binmdcint.o

#R4DIVO = four_caspt2_module.o nbitsa.o readvec.o read1mo.o \
	readorb_enesym.o \
	one_e_exct.o dim.o diag.o cutoff.o \
	readint2_ivo.o fockivo.o \
	calce0.o \
	takekr.o \
	intmo.o \
	timing.o mem.o \
	uramda_s_half.o nrintread.o \
	checkdgc.o e0test_v2.o casci.o casdet.o casmat.o r4divo.o

#R4DIVO_TY = four_caspt2_module.o nbitsa.o \
	read1mo_ty.o \
	readorb_enesym_ty.o \
	readint2_ivo_ty.o \
	fockivo_ty.o \
	intmo_ty.o \
	pgsym_ty.o \
	takekr.o dim.o one_e_exct.o diag.o \
	timing.o mem.o \
	uramda_s_half.o \
	e0test_v2.o r4divo_co.o
#	 cutoff.o readvec.o

#R4DCASPT2O = four_caspt2_module.o nbitsa.o readvec.o read1mo.o \
	readorb_enesym.o \
	one_e_exct.o dim.o diag.o cutoff.o \
	readint2_ord.o prtoutfock.o fockhf.o fockcasci.o \
	tramo.o e0after_tra.o trac.o sdet.o \
	calce0.o fockdiag.o \
	timing.o \
	intmo.o mem.o \
	intra.o \
	takekr.o \
	solvall_A_ord.o \
	solvall_B_ord.o \
	solvall_C_ord.o \
	solvall_D_ord.o \
	solvall_E_ord.o \
	solvall_F_ord.o \
	solvall_G_ord.o \
	solvall_H_ord.o \
	uramda_s_half.o nrintread.o \
	checkdgc.o e0test_v2.o casci.o casdet.o casmat.o r4dcaspt2_tra.o

#R4DCASPT2O_TY = four_caspt2_module.o nbitsa.o readvec.o \
	readorb_enesym_ty.o read1mo_ty.o \
	readint2_ord_ty.o \
	tramo_ty.o \
	one_e_exct.o dim.o diag.o cutoff.o \
	calce0.o \
	timing.o \
	mem.o \
	intra.o \
	takekr.o \
	solvall_A_ord_ty.o \
	solvall_B_ord_ty.o \
	solvall_C_ord_ty.o \
	solvall_D_ord_ty.o \
	solvall_E_ord_ty.o \
	solvall_F_ord_ty.o \
	solvall_G_ord_ty.o \
	solvall_H_ord_ty.o \
	uramda_s_half.o \
	pgsym_ty.o \
	checkdgc.o e0test_v2.o r4dcaspt2_tra_ty.o

R4DCASPT2O_CO = four_caspt2_module.o nbitsa.o readvec.o \
	readorb_enesym_co.o read1mo_co.o \
	readint2_ord_co.o \
	tramo_ty.o \
	one_e_exct.o dim.o diag.o cutoff.o \
	calce0.o \
	timing.o \
	mem.o \
	intra.o \
	takekr.o \
	solvall_A_ord_ty.o \
	solvall_B_ord_ty.o \
	solvall_C_ord_ty.o \
	solvall_D_ord_ty.o \
	solvall_E_ord_ty.o \
	solvall_F_ord_ty.o \
	solvall_G_ord_ty.o \
	solvall_H_ord_ty.o \
	uramda_s_half.o \
	checkdgc.o e0test_v2.o r4dcaspt2_tra_co.o

R4DIVO_CO = four_caspt2_module.o nbitsa.o \
        read1mo_co.o \
        readorb_enesym_co.o \
        readint2_ivo_ty.o \
        fockivo_ty.o \
        intmo_ty.o \
        pgsym_co.o \
        takekr.o dim.o one_e_exct.o diag.o \
        timing.o mem.o \
        uramda_s_half.o \
        e0test_v2.o r4divo_co.o
#        cutoff.o readvec.o

HFC_CASCI = four_caspt2_module.o nbitsa.o \
	one_e_exct.o dim.o hfc_casci.o

EEFF_CASCI = four_caspt2_module.o nbitsa.o \
	one_e_exct.o dim.o eeff_casci.o

NRMOBJS2 = four_caspt2_module.o nbitsa.o \
        readorb1_nr.o one_e_exct.o dim.o diag.o cutoff.o \
	readint2_nr.o prtoutfock.o fockhf.o fockcasci.o \
	tramo.o e0after_tra.o trac.o sdet.o \
	calce0.o fockdiag.o \
	solvall_A.o \
	solvall_B.o \
	solvall_C.o \
	solvall_D.o \
	solvall_E.o \
	solvall_F.o \
	solvall_G.o \
	solvall_H.o \
	uramda_s_half.o nrintread.o \
	checkdgc.o e0test_v2.o casci.o casdet.o casmat.o r4dcaspt2_ver2_nr.o

# This is a Intel mkl setting for the Institute for Molecular Science's linux server
# MKLROOT = /local/apl/lx/intel2020update2/compilers_and_libraries_2020.2.254/linux/mkl
MKLROOT = /local/apl/lx/intel2020update2/mkl
#BLASMOD = /local/apli/lx/intel2020update2/mkl/include/intel64/ilp64/blas95.mod
#LAPACKMOD = /local/apli/lx/intel2020update2/mkl/include/intel64/ilp64/lapack95.mod
#INC = -I$(BLASMOD) -I$(LAPACKMOD)
INC = -I$(MKLROOT)/include/intel64/ilp64 -i8 -I$(MKLROOT)/include
F90C = ifort
# F90FLAGS = $(INC) -mkl -DHAVE_ERF -FR -pad -O2 -mp1 -integer_size 64 -unroll
F90FLAGS = -mkl -DHAVE_ERF -FR -pad -O2 -mp1 -integer_size 64 -unroll


#all : r4divotyexe r4dcascityexe r4dcaspt2otyexe r4dcasciexe r4dcaspt2oexe r4divoexe hfc_casciexe eeff_casciexe
#all : r4dcasciexe r4dcaspt2oexe r4divoexe
all : r4divocoexe r4dcascicoexe r4dcaspt2ocoexe hfc_casciexe eeff_casciexe

#.f.o:
#	$(FORTRAN) $(OPTS) -c $*.f

.SUFFIXES: .f .o
.f.o:
	$(F90C) $(F90FLAGS) -c $*.f
#	$(F90C) $(F90FLAGS) -I$(MKLROOT)/include -c $*.f
#	$(F90C) $(F90FLAGS) -c $< :


#r4divotyexe : $(R4DIVO_TY)
#	$(F90C) $(F90FLAGS) -o $@ $(R4DIVO_TY) $(LAPACKLIB) $(BLASLIB)
#	mv r4divotyexe bin/r4divotyexe

#r4dcascityexe : $(R4DCASCI_TY)
#	$(F90C) $(F90FLAGS) -o $@ $(R4DCASCI_TY) $(LAPACKLIB) $(BLASLIB)
#	mv r4dcascityexe bin/r4dcascityexe

#r4dcasciexe : $(R4DCASCI)
#	$(F90C) $(F90FLAGS) -o $@ $(R4DCASCI) $(LAPACKLIB) $(BLASLIB)
#	mv r4dcasciexe bin/r4dcasciexe

#r4divoexe : $(R4DIVO)
#	$(F90C) $(F90FLAGS) -o $@ $(R4DIVO) $(LAPACKLIB) $(BLASLIB)
#	mv r4divoexe bin/r4divoexe

#r4dcaspt2oexe : $(R4DCASPT2O)
#	$(F90C) $(F90FLAGS) -o $@ $(R4DCASPT2O) $(LAPACKLIB) $(BLASLIB)
#	mv r4dcaspt2oexe bin/r4dcaspt2oexe

#r4dcaspt2otyexe : $(R4DCASPT2O_TY)
#	$(F90C) $(F90FLAGS) -o $@ $(R4DCASPT2O_TY) $(LAPACKLIB) $(BLASLIB)
#	mv r4dcaspt2otyexe bin/r4dcaspt2otyexe

hfc_casciexe : $(HFC_CASCI)
	$(F90C) $(F90FLAGS) -o $@ $(HFC_CASCI)
	mv hfc_casciexe bin/hfc_casciexe

eeff_casciexe : $(EEFF_CASCI)
	$(F90C) $(F90FLAGS) -o $@ $(EEFF_CASCI)
	mv eeff_casciexe bin/eeff_casciexe

r4dcascicoexe : $(R4DCASCI_CO)
	$(F90C) $(F90FLAGS) -o $@ $(R4DCASCI_CO)
	mv r4dcascicoexe bin/r4dcascicoexe

r4dcaspt2ocoexe : $(R4DCASPT2O_CO)
	$(F90C) $(F90FLAGS) -o $@ $(R4DCASPT2O_CO)
	mv r4dcaspt2ocoexe bin/r4dcaspt2ocoexe

r4divocoexe : $(R4DIVO_CO)
	$(F90C) $(F90FLAGS) -o $@ $(R4DIVO_CO)
	mv r4divocoexe bin/r4divocoexe

#nrmain2 : $(NRMOBJS2)
#	 $(FC) $(FFLAGS) -o $@ $(NRMOBJS2) $(LAPACKLIB) $(BLASLIB)
#

clean:
	rm *.o
	rm *.mod
