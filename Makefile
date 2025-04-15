OS := $(shell uname -s)

ifeq ($(OS), Windows_NT)
    CC = cl.exe
    CFLAGS = /EHsc /std:c++17 /openmp
    LDFLAGS = 
    RM = del /Q
else
    CC = g++
    CFLAGS = -std=c++17 -Wall -fopenmp
    LDFLAGS = -fopenmp
    RM = rm -f
endif

TARGET = Particle

SRC = Particle.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	$(RM) $(TARGET)

