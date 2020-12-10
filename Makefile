CC       = clang
CFLAGS   = -Wall -g
LDFLAGS  = -lm
OBJFILES = renderer.o obj_parser.o list.o
TARGET   = renderer

all: $(TARGET)

renderer: $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o renderer $(OBJFILES)

parser: obj_parser.o list.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o parser obj_parser.o list.o


renderer.o: renderer.c obj_parser.c list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c renderer.c obj_parser.c list.c

obj_parser.o: obj_parser.c list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c obj_parser.c list.c

list.o: list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c list.c

clean:
	rm -f $(OBJFILES) $(TARGET)
