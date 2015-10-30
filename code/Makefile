CC=g++
CFLAGS=-Werror
EXEC=travellingsp
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

travellingsp: $(OBJ)
	    @$(CC) -o $@ $^

%.o: %.c
	@$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	@rm -rf *.o

mrproper: clean
	@rm -rf $(EXEC)