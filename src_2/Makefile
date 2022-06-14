#gfortran -W MLPGR_frame1.f strmain.f90 regrid_old.for -o a

IFORTCHK=-check all -fpe0 -warn -traceback

COMPILER = ifort

ifeq ($(COMPILER), ifort)
	MODPATH=-module subroutines
	C_COMP=icc
	#PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64
	PARA_FLAGS=-lstdc++ -fopenmp -L/usr/local/cuda/lib64
	PARA_PATH=$(HOME)/Programs/paralution_icc_min
	OTHER_FLAGS = -integer-size 32
endif

PARA_CONV=$(PARA_PATH)/src/plug-ins/paralution_fortran.cpp
PARA_LIB=$(PARA_PATH)/build/lib/libparalution.so
#PARA_LIB=./libparalution.so
PARA_INC=-I$(PARA_PATH)/build/inc

MAIN=mlpgrMain
HEADER=mlpgrHeader

SUBM=\
subroutines/mlpgMainSubs.f

SUBS=\
subroutines/$(HEADER).f90 \
subroutines/modCommon.f90 \
subroutines/modules_v3.1.f90 \
subroutines/fnptCoupling.f90 \
subroutines/nodelinkNew_v2.3.f90 \
subroutines/interpNew_v1.5.f90 \
subroutines/interpFnc_v1.5.f90 \
subroutines/collision_v1.f90 \
subroutines/remesh_v1.f90 \
subroutines/resume.f90

OBJ=paralution_fortran.o $(HEADER).o $(SUBM) $(MAIN).f90

default: $(OBJ)
	mkdir -p Export Output OutputPara
	$(COMPILER) -o mlpgr $^ $(PARA_LIB) $(PARA_FLAGS) $(MODPATH) $(OTHER_FLAGS)

$(HEADER).o: $(SUBS)
	$(COMPILER) -c $< $(PARA_FLAGS) $(OTHER_FLAGS) $(MODPATH)

paralution_fortran.o: $(PARA_CONV)
	$(C_COMP) $(PARA_FLAGS) $(PARA_INC) -c $(PARA_CONV)

clean:
	rm -rf Export Output OutputPara
	mkdir -p Export Output OutputPara

cleanAll:
	rm -rf ./mlpgr ./outpara Export Output OutputPara subroutines/*.mod *.o

run: clean default
