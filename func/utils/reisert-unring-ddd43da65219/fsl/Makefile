include ${FSLCONFDIR}/default.mk

PROJNAME=unring

USRINCFLAGS=-I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS=-L${LIB_NEWMAT} 

OBJS=unring.o
LIBS=-lwarpfns -lmeshclass -lbasisfield -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lutils -lnewmat -lm -lz -lfftw3

XFILES=unring 
FXFILES=

all: ${XFILES} schedule

schedule:
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi

unring: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}
