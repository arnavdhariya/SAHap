CC=gcc -std=c99 -O2 -Ilibwayne/include -DNDEBUG=1 -Wno-unused-variable -Wno-unused-but-set-variable -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow

sahap-c: sahap-c.o
	$(CC) -o sahap-c *.o -lm -Llibwayne -lwayne -lm

clean:
	rm -f *.o *.exe*
