FC	= mpif90
LD	= $(FC)
RM	= /bin/rm -f
OLEVEL	= -O3
DBG := -g -fbounds-check -Wall -fcheck=all
FOPTS	= -mcmodel=medium -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -fbounds-check -Wall -cpp  #-std=f2008ts -g
LIB = -llapack   							#for deigo
#FOPTS	= -align none -mcmodel medium -warn all
FFLAGS	= $(FOPTS) $(OLEVEL) $(DBG)

FLWOBJS = \
./common_data.f90 \
./flood_fill.f90 \
./paraview.f90 \
./get_interface.f90\
./read_fields.f90 \
./main.f90 

OBJS	= $(FLWOBJS)
EXEC    =  ./drop_count

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) -o $@ $(LIB)

clean:
	$(RM) $(EXEC)
	rm -f *.mod

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
