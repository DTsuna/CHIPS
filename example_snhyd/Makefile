FC	      = gfortran

FFLAGS        = -fno-automatic

DEST	      = ${HOME}/bin

EXTHDRS	      =

HDRS	      =

LDFLAGS	      =

LIBS	      =

LINKER	      = gfortran

MAKEFILE      = Makefile

EXT_LIST      = .f .f90

EXT           = .f

OBJDIR        =	obj

NO            = src/inclcnst.f src/inclion.f src/inclmn.f src/inclm1.f src/inclold.f src/view.f src/opac_original.f src/eos_ionize.f

CPPS          =	$(shell find * -name '*$(EXT)')
SRCS          =	$(filter-out $(NOMAKE) $(NOFILE), $(CPPS))
DIRS          =	$(dir $(SRCS))
BINDIRS       =	$(addprefix $(OBJDIR)/, $(DIRS))
NOOBJ         = $(addprefix $(OBJDIR)/, $(patsubst %.f, %.o, $(NO)))

OBJS          =	$(filter-out $(NOOBJ), $(addprefix $(OBJDIR)/, $(patsubst %$(EXT), %.o, $(SRCS))))
OBJA          = $(shell find * -name '*.o')


PRINT	      = pr

PROGRAM	      = runsnhyd



default:
		@[ -d  $(OBJDIR)   ] || mkdir -p $(OBJDIR)
		@[ -d "$(BINDIRS)" ] || mkdir -p $(BINDIRS)
		@make $(OBJS) --no-print-directory
	

all:		$(PROGRAM)

$(PROGRAM):     $(OBJA)
		$(LINKER) $(LDFLAGS) $(OBJA) $(LIBS) -o $(PROGRAM)

$(OBJDIR)/%.o: %$(EXT)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:;		rm -rf $(OBJS) $(PROGRAM) $(OBJDIR)

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

print:;		$(PRINT) $(HDRS) $(SRCS)

program:        $(PROGRAM)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@make -f $(MAKEFILE) DEST=$(DEST) install
###
advanc.o: inclm1.f inclmn.f inclold.f
alar.o: inclm1.f inclmn.f
atmos.o: inclm1.f inclmn.f
avrge.o: inclm1.f inclmn.f
cournt.o: inclm1.f inclmn.f
deltaa.o: inclm1.f inclmn.f
eos.o: inclm1.f inclmn.f inclcnst.f inclion.f
extend.o: inclm1.f inclmn.f
grav.o: inclm1.f inclmn.f
grow.o: inclm1.f inclmn.f
init.o: inclm1.f inclmn.f
intfce.o: inclm1.f inclmn.f
old.o: inclm1.f inclmn.f inclold.f
opac.o: inclm1.f inclmn.f inclcnst.f
output.o: inclm1.f inclmn.f
radtra.o: inclm1.f inclmn.f inclcnst.f inclold.f
riemnt.o: inclm1.f inclmn.f
saha.o: inclm1.f inclmn.f inclcnst.f inclion.f
snhyd.o: inclm1.f inclmn.f
state.o: inclm1.f inclmn.f
view.o: inclmn.f
