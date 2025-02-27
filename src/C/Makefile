CC=gcc -O3 #-ggdb

sahap-c: bitvec.o sahap-c.o misc.o sets.o
	$(CC) -o sahap-c *.o -lm -Llibwayne -lwayne

clean:
	rm -f *.o *.exe*
