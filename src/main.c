#include <stdio.h>
#include <stdlib.h>
#include "hexdump.h"

#define ESC 0x1b

int main(int argc, char **argv)
{
	FILE *prn;

	// check command line args
	if (argc < 2) {
		printf("ERROR: must specify device name\n");
		return -1;
	}

	// open printer device
	if ((prn = fopen(argv[1], "r+b")) == NULL) {
		printf("ERROR: couldn't open printer device '%s'\n", argv[1]);
		return -1;
	}

	// INITIALISE
	fprintf(prn, "%c%c", ESC, '@');

	// REQUEST STATUS
	fprintf(prn, "%c%c%c", ESC, 'i', 'S');

	// Read status buffer from printer
	unsigned char buf[32];
	int timeout = 128;
	do {
		fread(buf, 1, 32, prn);
	} while ((buf[0] != 0x80) && (timeout-- > 0));

	if (timeout > 0) {
		printf("Printer status:\n");
		hex_dump(buf, 32);
	} else {
		printf("TIMEOUT\n");
		return -1;
	}

	// Close the printer stream
	fclose(prn);
	return 0;
}
