# https://stackoverflow.com/questions/30573481/path-include-and-src-directory-makefile

TARGET   = renderer
OBJ_DIR  = obj
SRC_DIR  = src

SRC     := $(wildcard $(SRC_DIR)/*.c)
OBJ 	:= $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC))

CC       = clang
CFLAGS   = -Wall -g -Iinclude
OPTFLAGS = -O3 -Ofast -msse3 -mfpmath=sse
LDFLAGS  = -lm -lGL -lGLU -lglut -lpthread
LDLIBS   =

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

clean:
	$(RM) -rv $(OBJ_DIR) $(TARGET)

-include $(OBJ:.o=.d)
