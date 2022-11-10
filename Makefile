CC = c++
CFLAGS = -Wall -g
FLAGS = -Wall -g
LIBS = -lm
OBJS = SCP.o mmas_ml_main.o pheromone.o


mmas_ml: $(OBJS)
	$(CC) $(FLAGS) -o mmas_ml $(OBJS) $(LIBS)
.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	/bin/rm -rf *.o *~ mmas_ml $(OBJS) $(TARGET)
