BINS = asmuthbloom

CFLAGS = -std=c99 -g -Wall -Wextra -Wno-unused -Wno-sign-compare

.PHONY: all clean

all: $(BINS)

clean:
	$(RM) $(BINS)

asmuthbloom: LDLIBS = -lgmp
