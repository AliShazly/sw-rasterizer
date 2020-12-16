CC       = clang
CFLAGS   = -Wall -g
LDFLAGS  = -lm -lGL -lGLU -lglut
OBJFILES = display.o renderer.o obj_parser.o list.o wu_line.o
TARGET   = renderer

all: $(TARGET)

renderer: $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o renderer $(OBJFILES)

display.o: display.c renderer.c obj_parser.c list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c display.c renderer.c obj_parser.c list.c

renderer.o: renderer.c obj_parser.c list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c renderer.c obj_parser.c list.c

obj_parser.o: obj_parser.c list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c obj_parser.c list.c

wu_line.o: wu_line.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c wu_line.c

list.o: list.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c list.c

clean:
	rm -f $(OBJFILES) $(TARGET)
