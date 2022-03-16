CC = g++

CFLAGS = -g
LDFLAGS = 

REQ = vec5d.o eostable.o RK45.o main.o

%.o: %.cpp %.hpp    ; $(CC) $(CFLAGS) -c $<

main: $(REQ)    ; $(CC) $(CFLAGS) $(REQ) $(LDFLAGS) -o $@

clean:  ; rm -f *.o main
