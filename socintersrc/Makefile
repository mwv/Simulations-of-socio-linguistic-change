objects = socinterfuncs.o socinter.o

socinter : $(objects)
	gcc -o socinter -lm -O3 -Wall -Werror $(objects)
socinterfuncs.o : socinterfuncs.c
	gcc -c -Wall -Werror -O3 socinterfuncs.c
socinter.o : socinter.c
	gcc -c -Wall -Werror -O3 socinter.c
clean : 
	rm socinter $(objects)