FC	      = gfortran

FFLAGS        = -fno-automatic

DEST	      = ${HOME}/bin

EXTHDRS	      =

HDRS	      =

LDFLAGS	      =

LIBS	      =

LINKER	      = gfortran

MAKEFILE      = Makefile

OBJS	      = snhyd/advanc.o \
		snhyd/alar.o \
		snhyd/atmos.o \
		snhyd/avrge.o \
		snhyd/compo.o \
		snhyd/cournt.o \
		snhyd/deltaa.o \
		snhyd/eos.o \
		snhyd/extend.o \
		snhyd/grav.o \
		snhyd/grow.o \
		snhyd/init.o \
		snhyd/intfce.o \
		snhyd/old.o \
		snhyd/opac.o \
		snhyd/output.o \
		snhyd/radtra.o \
		snhyd/riemnt.o \
		snhyd/saha.o \
		snhyd/snhyd.o \
		snhyd/state.o \
		snhyd/eos_helm.o \
		snhyd/eos_helm_e.o \
                snhyd/eos_helm_p.o \
		snhyd/tote.o

PRINT	      = pr

PROGRAM	      = runsnhyd

SRCS	      = snhyd/advanc.f \
		snhyd/alar.f \
		snhyd/atmos.f \
		snhyd/avrge.f \
		snhyd/compo.f \
		snhyd/cournt.f \
		snhyd/deltaa.f \
		snhyd/eos.f \
		snhyd/extend.f \
		snhyd/grav.f \
		snhyd/grow.f \
		snhyd/inclcnst.f \
		snhyd/inclion.f \
		snhyd/inclm1.f \
		snhyd/inclmn.f \
		snhyd/inclold.f \
		snhyd/init.f \
		snhyd/intfce.f \
		snhyd/old.f \
		snhyd/opac.f \
		snhyd/output.f \
		snhyd/radtra.f \
		snhyd/riemnt.f \
		snhyd/saha.f \
		snhyd/snhyd.f \
		snhyd/state.f \
		snhyd/eos_helm.f90 \
		snhyd/eos_helm_e.f90 \
                snhyd/eos_helm_p.f90 \
		snhyd/tote.f

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f $(OBJS)

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
