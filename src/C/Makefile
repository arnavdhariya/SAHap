CC=gcc -std=c99 -O2 -DNDEBUG=1 -Wno-unused-variable -Wno-unused-but-set-variable -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow

sahap-c: bitvec.o sahap-c.o misc.o sets.o
	$(CC) -o sahap-c *.o -lm -Llibwayne -lwayne

clean:
	rm -f *.o *.exe*
