BINS = asmuthbloom

CFLAGS = -g -Wall -Wextra

.PHONY: all clean

all: $(BINS)

clean:
	$(RM) $(BINS)

asmuthbloom: LDLIBS = -lgmp
