CC = gcc
CFLAGS = -U__STRICT_ANSI__ -Wall #-pthread #-pg -fno-inline

all: scripts/delete-commas-inside-quotes

scripts/delete-commas-inside-quotes: scripts/delete-commas-inside-quotes.c
	gcc -o scripts/delete-commas-inside-quotes scripts/delete-commas-inside-quotes.c

clean:
	rm -f scripts/delete-commas-inside-quotes
