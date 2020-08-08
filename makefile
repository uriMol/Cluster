FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main 

main: spmat.o main.o
	gcc main.o spmat.o -o prog $(LIBS)
spmat.o: spmat.c
	gcc $(FLAGS) -c spmat.c
main.o: main.c
	gcc $(FLAGS) -c main.c

clean:
	rm -rf *.o prog
