all: main
	@echo Para executar: ./main sistema1x3x3.dat ou outro .dat
main: main.o Algoritmos.o
	gcc main.o Algoritmos.o -o main

main.o: main.c
	gcc -c main.c

Algoritmos.o: Algoritmos.c
	gcc -c Algoritmos.c

clean:
	$(RM) main$(EXT) main.o Algoritmos.o
	clear