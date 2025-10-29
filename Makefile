FUNC := clang++
copt := -c
OBJ_DIR := ./bin/
FLAGS := -O3 -lm -g -Werror -march=native -flto -ftree-vectorize -ffast-math

C_FILES := $(wildcard src/*.cc)
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(C_FILES:.c=.obj)))

all:
	$(FUNC) ./main.cc -o ./main.exe $(FLAGS)

clean:
	rm -f ./*.exe