FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm

all: main 

main: spmat.o spmatUtils.o list.o main.o group.o modularityMax.o mainUtils.o
	gcc main.o spmat.o spmatUtils.o list.o group.o modularityMax.o mainUtils.o -o prog $(LIBS)
mainUtils.o: mainUtils.c
	gcc $(FLAGS) -c mainUtils.c
list.o:	list.c
	gcc $(FLAGS) -c list.c
spmatUtils.o:	spmatUtils.c
	gcc $(FLAGS) -c spmatUtils.c
spmat.o: spmat.c
	gcc $(FLAGS) -c spmat.c
modularityMax.o: modularityMax.c
	gcc $(FLAGS) -c modularityMax.c
group.o: group.c
	gcc $(FLAGS) -c group.c
main.o: main.c
	gcc $(FLAGS) -c main.c
clean:
	rm -rf *.o prog
