CC = c++
CFLAGS = -Wall
FLAGS = -Wall -O2
LIBS = -lm
OBJS = SCPv.o mmas_ml_main.o


mmas_ml: $(OBJS)
	$(CC) $(FLAGS) -o mmas_ml $(OBJS) $(LIBS)
.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	/bin/rm -rf *.o *~ mmas_ml $(OBJS) $(TARGET)
