CC = gcc

dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

HTSLIB_DIR = ../htslib

ifeq ($(DEBUG),1)
    CFLAGS = -O0 -g -Wall -Wextra -DDEBUG -Iinclude -Iklib -I. -I$(HTSLIB_DIR) -fopenmp
else
    CFLAGS = -O3 -Wall -Wextra -Iinclude -Iklib -I. -I$(HTSLIB_DIR) -fopenmp
endif

LDFLAGS = -L$(HTSLIB_DIR) -fopenmp \
  -Wl,-rpath,'$$ORIGIN/../../htslib'
LDLIBS = -lm -lhts -lz

SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c,obj/%.o,$(SRC))

BIN = bin/GLIMPSE2_build_ser_from_bcf

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS)

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f obj/*.o $(BIN)

.PHONY: all clean
