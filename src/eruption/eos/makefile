FC	:=	gfortran
FFLAGS	:=	-fno-automatic -fno-range-check


#You should not change the following
OBJDIR	:=	obj

CPPS	:=	$(shell find * -name '*.f90')
SRCS	:=	$(filter-out $(NOMAKE), $(CPPS))
DIRS	:=	$(dir $(SRCS))
BINDIRS	:=	$(addprefix $(OBJDIR)/, $(DIRS))

NOFILE	:= 	inclcnst.f inclion.f inclmn.f inclm1.f inclold.f view.f opac_original.f eos_ionize.f
NOOBJ	:= 	$(addprefix $(OBJDIR)/, $(patsubst %.f90, %.o, $(NOFILE)))

OBJS	:=	$(filter-out $(NOOBJ), $(addprefix $(OBJDIR)/, $(patsubst %.f90, %.o, $(SRCS))))

ifeq "$(strip $(OBJDIR))" ""
	OBJDIR = .
endif

ifeq "$(strip $(DIRS))" ""
	OBJDIR = .
endif

default:
	@[ -d  $(OBJDIR)   ] || mkdir -p $(OBJDIR)
	@[ -d "$(BINDIRS)" ] || mkdir -p $(BINDIRS)
	@make $(OBJS) --no-print-directory

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) $(DEPS)

remake:
	make clean && make

-include $(DEPS)
