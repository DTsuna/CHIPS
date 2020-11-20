FC	:=	gfortran
FFLAGS	:=	-fno-automatic -fno-range-check
PROG	:=	runsnhyd

CPPS	=	$(shell find * -name '*.o')

OBJS	=	$(CPPS)


all:	$(OBJS) $(PROG)
	$(eval OBJS = $(shell find * -name '*.o'))
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(PROG)

$(PROG):
	@cd f; make; cd ../f90; make; cd ..

clean:
	rm -f $(PROG) $(OBJS)
