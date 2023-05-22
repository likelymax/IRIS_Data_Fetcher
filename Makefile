find_event: find_event.c find_event.o
	gcc -o find_event find_event.o

find_event.o: find_event.c
	gcc -c find_event.c
	
clean:
	rm -f *.o
