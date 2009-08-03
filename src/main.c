/**************************
 * P-Touch PT-2450DX printer driver
 *
 * P. Pemberton, 2009
 *
 * Specs:
 *   Printer head is 128 dots, 180dpi, for a total print area of ~18mm vertical
 *   by however long your TZ label tape is.
 *   Printhead size is (128/180)=0.711[.] inches, or 18.0622[.] mm
 *   Each dot is (18.062/128) = 0.1411[.] mm
 *   Printable area is (numdots * 0.1411) mm
 *
 *     Tape width   Margins      Printable area
 *     mm    dots   mm     dots  mm    dots
 *      6mm   42    1.0mm     7   4mm   28
 *      9mm   63    1.0mm     7   7mm   49
 *     12mm   85    2.0mm    14   8mm   57
 *     18mm  127    3.0mm    21  12mm   85
 *     24mm  170    3.0mm   128  18mm  128 ***
 *
 *   24mm is slightly odd. Because the printhead is only 128 dots (18mm), the
 *   margins are enforced by this, and not the driver software. It is impossible
 *   to print right to the edge of a 24mm label in a PT-2450DX.
 *
 **************************/

#include <stdio.h>
#include <stdlib.h>
#include "ptouch.h"

/****************************************************************************/

int main(int argc, char **argv)
{
	pt_Device *dev;

	// check command line args
	if (argc < 2) {
		printf("ERROR: must specify device name\n");
		return -1;
	}

	// Open and initialise the printer
	dev = pt_Initialise(argv[1]);

	if (dev == NULL) {
		printf("Error opening printer device.\n");
		return -1;
	}

	// Close the printer device
	pt_Close(dev);

	return 0;
}
