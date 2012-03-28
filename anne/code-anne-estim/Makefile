SRCS = alloc.c logL.c prior.c moves.c mcmc.c init.c InputOutput.c tuneVariances.c main.c
OBJS = alloc.o logL.o prior.o moves.o mcmc.o init.o InputOutput.o tuneVariances.o main.o
EXE = estimation

.c.o: 
	gcc -g3 -Wall -c -fmessage-length=0 $<

all:	$(OBJS)
	gcc -o ${EXE} -O3 ${OBJS} -lgsl -lgslcblas -lm --static

clean:
	rm *.exe *.o
