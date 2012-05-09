SRCS = matvec.c alloc.c output.c main.c
OBJS = matvec.o alloc.o output.o main.o
EXE = simulation

.c.o: 
	gcc -g3 -Wall -c -fmessage-length=0 $<

all:	$(OBJS)
	gcc -o ${EXE} -O3 ${OBJS} -lgsl -lgslcblas -lm --static

clean:
	rm *.exe *.o
