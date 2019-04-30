CC = gcc
CFLAGS = -fPIC -O3

DEPS = spherical.h
OBJ  = spherical.o

default: libspherical.a

libspherical.a:
	$(CC) $(CFLAGS) -c spherical.c
	ar rcs $@ $(OBJ)

clean:
	rm *.o *.a
