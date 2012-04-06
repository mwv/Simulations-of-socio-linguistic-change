objects = socimpactfuncs.o socimpact.o

socimpact : $(objects)
	gcc -o socimpact -lm -O3 -Wall -Werror $(objects)
socimpactfuncs.o : socimpactfuncs.c
	gcc -c -O3 -Wall -Werror socimpactfuncs.c
socimpact.o : socimpact.c
	gcc -c -O3 -Wall -Werror socimpact.c
clean :
	rm socimpact $(objects)