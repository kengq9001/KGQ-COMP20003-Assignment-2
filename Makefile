voronoi2: main.o wt_ops.o dcel_ops.o output.o
	gcc -Wall -o voronoi2 main.o wt_ops.o dcel_ops.o output.o -lm -g

main.o: main.c wt_ops.h dcel_ops.h output.h
	gcc -Wall -o main.o main.c -c -g

wt_ops.o: wt_ops.c wt_ops.h
	gcc -Wall -o wt_ops.o wt_ops.c -c -g

dcel_ops.o: dcel_ops.c dcel_ops.h
	gcc -Wall -o dcel_ops.o dcel_ops.c -c -g

output.o: output.c output.h wt_ops.h dcel_ops.h
	gcc -Wall -o output.o output.c -c -g

clean: voronoi2
	rm *.o voronoi2
