PROG	:=	a.out
CC	:=	icc
CFLAGS	:=	-Wall -ip -ipo -parallel -qopenmp -O3# -fbounds-check -fstack-check

INC	:=	-I./include
NOMAKE	:=	outp-data% inp-data% nomakedir% transfer%

#You should not change the following
OBJDIR	:=	obj

CPPS	:=	$(shell find * -name '*.c')
SRCS	:=	$(filter-out $(NOMAKE), $(CPPS))
DIRS	:=	$(dir $(SRCS))
BINDIRS	:=	$(addprefix $(OBJDIR)/, $(DIRS))

OBJS	:=	$(addprefix $(OBJDIR)/, $(patsubst %.c, %.o, $(SRCS)))
DEPS	:=	$(OBJS:%.o=%.d)
TILS	:=	$(patsubst %.c, %.c~, $(SRCS))

ifeq "$(strip $(OBJDIR))" ""
	OBJDIR = .
endif

ifeq "$(strip $(DIRS))" ""
	OBJDIR = .
endif

default:
	@[ -d  $(OBJDIR)   ] || mkdir -p $(OBJDIR)
	@[ -d "$(BINDIRS)" ] || mkdir -p $(BINDIRS)
	@make object --no-print-directory

all:	$(OBJS) $(PROG)

object:	$(OBJS)

$(PROG):	$(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(INC) -o $@ $^ -lm

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INC) -o $@ -c -MMD -MP $<

clean:
	rm -rf $(PROG) $(OBJDIR) $(DEPS)

remake:
	make clean && make

-include $(DEPS)
