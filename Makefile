CC       = clang
CFLAGS   = -Wall -g
LDFLAGS  = -lm -lGL -lGLU -lglut -lpthread
OBJFILES = display.o renderer.o obj_parser.o list.o wu_line.o
TARGET   = renderer

all: $(TARGET)

renderer: $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJFILES)

clean:
	rm -f $(OBJFILES) $(TARGET)
