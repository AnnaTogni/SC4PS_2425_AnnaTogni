CC     = gcc
CFLAGS = -Wall -Wextra -O2
LDLIBS = -lm

numericalsum_v2: NumericalSum_v2.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f numericalsum_v2

.PHONY: clean
