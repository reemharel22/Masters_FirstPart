all: soNumerit soFL soEF so

soNumerit: SuOlsonNumerit.c  
	gcc -o soNumerit SuOlsonNumerit.c -lm -fopenmp

soEF: SuOlsonEddingtonFactor.c
	gcc -o soEF SuOlsonEddingtonFactor.c -lm -fopenmp

so: tridFunc.o SuOlsonAll.c
	gcc -o so SuOlsonAll.c tridFunc.o -lm -fopenmp

tridFunc: tridFunc.c tridFunc.h
	gcc -c -o tridFunc.o tridFunc.c -lm -fopenmp

clean:
	rm -f *.o soNumerit
