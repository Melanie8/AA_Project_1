CC=gcc
CFLAGS=-Werror
EXEC=travellingsp
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)

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