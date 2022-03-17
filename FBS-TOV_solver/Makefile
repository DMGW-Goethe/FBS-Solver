CC = g++

CFLAGS = -g -O3
LDFLAGS =

REQ = vector.o eos.o RK45.o main.o

%.o: %.cpp %.hpp    ; $(CC) $(CFLAGS) -c $<

main: $(REQ)    ; $(CC) $(CFLAGS) $(REQ) $(LDFLAGS) -o $@

clean:  ; rm -f *.o main
