FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main 

main: spmat.o main.o
	gcc main.o spmat.o list.o -o prog $(LIBS)
list.o:	list.c
	gcc $(FLAGS) -c list.c
spmat.o: spmat.c
	gcc $(FLAGS) -c spmat.c
main.o: main.c
	gcc $(FLAGS) -c main.c


clean:
	rm -rf *.o prog
