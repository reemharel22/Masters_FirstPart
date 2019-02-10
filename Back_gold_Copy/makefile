all: so

so: tridFunc.o SuOlsonAll.c SuOlsonAll.h
	gcc -o so SuOlsonAll.c tridFunc.o -lm -fopenmp

tridFunc: tridFunc.c tridFunc.h
	gcc -c -o tridFunc.o tridFunc.c -lm -fopenmp

clean:
	rm -f *.o soNumerit
