CC     = gcc
CFLAGS = -Wall -Wextra -O2
LDLIBS = -lm

numericalsum: NumericalSum.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f numericalsum

.PHONY: clean
