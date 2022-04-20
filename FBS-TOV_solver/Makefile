CC = g++

CFLAGS = -g -O3 -I/usr/include/python3.8
LDFLAGS = -lpython3.8

REQ = vector.o eos.o RK45.o nsmodel.o plotting.o

all: $(REQ)
	$(CC) $(CFLAGS) $(REQ) main.cpp -o main $(LDFLAGS)


%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $<


clean:
	rm -f *.o main
