FC	:=	gfortran
FFLAGS	:=	-fno-automatic -fno-range-check
PROG	:=	eruption

TARGETS	:=	src/eruption/hydro src/eruption/eos	

CPPS	=	$(shell find $(TARGETS) -name '*.o')

OBJS	=	$(CPPS)


all:	$(OBJS) $(PROG)
	$(eval OBJS = $(shell find $(TARGETS) -name '*.o'))
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $(PROG)

$(PROG):
	@for target in $(TARGETS); do \
	(\
		cd $$target; \
		make; \
	)\
	done
	cd src/LC; \
	python3 lcsetup.py install --user

clean:
	rm -f $(PROG) $(OBJS)
	@cd src/LC; \
	python3 lcsetup.py clean --all
