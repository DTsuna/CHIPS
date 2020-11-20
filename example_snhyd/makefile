FC	:=	gfortran
FFLAGS	:=	-fno-automatic -fno-range-check
PROG	:=	runsnhyd

TARGETS	:=	f f90

CPPS	=	$(shell find * -name '*.o')

OBJS	=	$(CPPS)


all:	$(OBJS) $(PROG)
	$(eval OBJS = $(shell find * -name '*.o'))
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(PROG)

$(PROG):
	@for target in $(TARGETS); do \
	(\
		cd $$target; \
		make; \
	)\
	done

clean:
	rm -f $(PROG) $(OBJS)
