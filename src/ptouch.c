/****************************************************************************
 * Project:   P-touch printer driver library
 * Developer: Philip Pemberton
 * Purpose:   Make Brother P-touch (PT-series) printers do something besides
 *            gather dust.
 *
 *            Currently supports:
 *              PT-2450DX
 ****************************************************************************/

// TODO: disable
#define DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hexdump.h"
#include "pt_image.h"
#include "ptouch.h"

#define ESC 0x1b

pt_Device *pt_Initialise(char *path)
{
	pt_Device	*dev;
	FILE		*prn;

	// Try and open the printer device
	prn = fopen(path, "r+b");
	if (prn == NULL) {
		return NULL;
	}

	// Printer device open, send an init command and read the status
	fprintf(prn, "%c%c", ESC, '@');

	// Allocate memory for the device block
	dev = malloc(sizeof(pt_Device));
	if (dev == NULL) {
		fclose(prn);
		return NULL;
	}

	// Store the file pointer
	dev->fp = prn;

	// Memory allocation OK, now get the printer's status
	if (pt_GetStatus(dev) == 0) {
		return dev;
	} else {
		free(dev);
		return NULL;
	}
}

int pt_GetStatus(pt_Device *dev)
{
	// REQUEST STATUS
	fprintf(dev->fp, "%c%c%c", ESC, 'i', 'S');

	// Read status buffer from printer
	unsigned char buf[32];
	int timeout = 128;
	do {
		fread(buf, 1, 32, dev->fp);
	} while ((buf[0] != 0x80) && (timeout-- > 0));

	// Check for timeout
	if (timeout == 0) {
		// Timeout
		return -1;
	}

#ifdef DEBUG
	printf("DEBUG: Printer status buffer = \n");
	hex_dump(buf, 32);
#endif

	// Decode the status buffer, store the results in the device object
	dev->headMark     = buf[0];
	dev->size         = buf[1];
	dev->errorInfo[0] = buf[8];
	dev->errorInfo[1] = buf[9];
	dev->mediaWidth   = buf[10];
	dev->mediaType    = buf[11];
	dev->mediaLength  = buf[17];
	dev->statusType   = buf[18];
	dev->phaseType    = buf[19];
	dev->phaseHi      = buf[20];
	dev->phaseLo      = buf[21];
	dev->notification = buf[22];

	// Operation succeeded
	return 0;
}

// TODO: print options struct parameter (e.g. fullcut, halfcut, print res,
//
int pt_Print(pt_Device *dev, pt_Image *image)
{
	// TODO: trap dev == NULL
	// TODO: trap image == NULL
	// TODO: trap image.height > printhead.height
	// TODO: trap image.height <= 0
	// TODO: trap image.width <= 0
	//
	// allocate print buffer
	//
	// pack pixels -- 8 pixels => 1 byte
	//
	// compress print data (packbits)
	//
	// send print buffer to printer
	//
	// free print buffer
}

void pt_Close(pt_Device *dev)
{
	// Sanity check -- make sure dev is not null
	if (dev == NULL) return;

	// Close the printer stream
	fclose(dev->fp);

	// Release the memory allocated to the status buffer
	free(dev);
}

