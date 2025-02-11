OS := $(shell uname -s)

ifeq ($(OS), Windows_NT)
    CC = cl.exe
    CFLAGS = /EHsc /std:c++17
    LDFLAGS = 
    RM = del /Q
else
    CC = g++
    CFLAGS = -std=c++17 -Wall
    LDFLAGS =
    RM = rm -f
endif

TARGET = Particle

SRC = Particle.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	$(RM) $(TARGET)
