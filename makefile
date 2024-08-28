CC = gcc

CFLAGS = \
	 -g \
	 -Iinclude \
	 -Wall \
	 -Werror \
	 -pedantic \
	 -fsanitize=address \
	 -fno-omit-frame-pointer \
     -fsanitize=float-cast-overflow \
	 -fopenmp 

CFLAGS_OPT = \
	     -O3 \
	     -march=native 

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o \
      obj/lattice.o \
      obj/potential.o \
      obj/run.o \
	  obj/verlet.o \
	  obj/fft.o \
	  obj/t1.o \
	  obj/t2.o \
	  obj/t3.o \

	  
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

H1: obj _H1

_H1: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ -o H1 $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj obj_test
	rm -rf *Zone.Identifier
	rm -f H1

.PHONY: clean
