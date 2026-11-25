CC ?= cc
CFLAGS ?= -O2 -fopenmp
all: hsweep hsweep2 hsweep_fast
%: %.c
	$(CC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f hsweep hsweep2 hsweep_fast
