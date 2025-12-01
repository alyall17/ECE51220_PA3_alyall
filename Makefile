CC = gcc
CFLAGS = -std=c99 -pedantic -Wvla -Wall -Wshadow -O3

all: pa3

pa3: pa3.c
	$(CC) $(CFLAGS) pa3.c -o pa3

clean:
	rm -f pa3
