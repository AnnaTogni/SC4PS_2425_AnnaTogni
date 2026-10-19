CC := gcc
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic -O2
GLIB_CFLAGS := $(shell pkg-config --cflags glib-2.0)
GLIB_LIBS := $(shell pkg-config --libs glib-2.0)

ifneq ($(wildcard /usr/include/hdf5/serial/hdf5.h),)
HDF5_CFLAGS ?= -I/usr/include/hdf5/serial
HDF5_LIBS ?= -lhdf5_serial
else ifneq ($(wildcard /usr/include/hdf5.h),)
HDF5_CFLAGS ?=
HDF5_LIBS ?= -lhdf5
else
HDF5_CFLAGS ?=
HDF5_LIBS ?=
endif

TARGET := homework6
SOURCES := main.c

.PHONY: all clean run run-config

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) $(GLIB_CFLAGS) $(HDF5_CFLAGS) -o $@ $(SOURCES) $(GLIB_LIBS) $(HDF5_LIBS) -lm

run: $(TARGET)
	./$(TARGET)

run-config: $(TARGET)
	./$(TARGET) --config homework6.cfg

clean:
	rm -f $(TARGET) output_chunks.h5
