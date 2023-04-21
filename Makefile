#make file - this is a comment section

CC=gcc  #compiler
CFLAGS=-I.
DEPS = g_power_x.h
OBJ = g_power_x.o functions.o
TARGET=output

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET):	$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
	

clean:
	rm *.o $(TARGET)

