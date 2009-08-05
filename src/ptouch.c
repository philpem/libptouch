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
//#define DEBUG

// This debug option forces Request Status to always "see" a good status
// block. Mostly useful for testing using write-to-file mode.
//#define DEBUG_SKIP_STATUS_READ

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gd.h>
#include "hexdump.h"
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
	if (pt_GetStatus(dev) != PT_ERR_SUCCESS) {
		// GetStatus failed, close the device and exit.
		pt_Close(dev);
		return NULL;
	}

	// Set printing parameters to defaults --
	//   Mirror off
	//   Autocut off
	dev->mirror = false;
	dev->autocut = false;

	return dev;
}

int pt_GetStatus(pt_Device *dev)
{
#ifdef DEBUG_SKIP_STATUS_READ
	unsigned char buf[32];
	memset(buf, 0x00, 32);
	buf[0] = 0x80;
	buf[1] = 0x20;
	buf[2] = 0x42;
	buf[3] = 0x30;
	buf[4] = 0x4b;
	buf[5] = 0x30;
	buf[10] = 0x0c;
	buf[11] = 0x01;
#else
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
		return PT_ERR_TIMEOUT;
	}
#endif

#ifdef DEBUG
	printf("DEBUG: Printer status buffer = \n");
	hex_dump(buf, 32);
#endif

	// Decode the status buffer, store the results in the device object
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

	// Set pixel width (label width in pixels)
	if (dev->mediaWidth >= 24) {
		// Label tape is 24mm or wider. Print head is 128 dots at 180dpi,
		// which is 18.06mm. Thus, we can only print on the centre 18mm
		// of a tape that is wider than 18mm.
		dev->pixelWidth = 128;
	} else {
		// Print head is 180dpi. Pixel size is mediaWidth * dpi. If we
		// multiply by ten, then divide by 254, we can avoid using
		// floating point to convert from inches to mm. The -2 is a
		// safety margin -- one pixel on either side of the label.
		// This is far closer than Brother suggest, but hey-ho.
		dev->pixelWidth = ((dev->mediaWidth * 180 * 10) / 254) - 2;
	}

	// Operation succeeded
	return PT_ERR_SUCCESS;
}

int pt_SetOption(pt_Device *dev, PT_E_OPTION option, int value)
{
	// trap dev == NULL
	if (dev == NULL) {
		return PT_ERR_BAD_PARAMETER;
	}

	// set option
	switch(option) {
		case PT_OPTION_MIRROR:			// Mirror
			dev->mirror = (value ? 1 : 0);
			return PT_ERR_SUCCESS;

		case PT_OPTION_AUTOCUT:		// Auto-cutter enable/disable
			dev->autocut = (value ? 1 : 0);
			return PT_ERR_SUCCESS;

		default:
			return PT_ERR_BAD_PARAMETER;
	}
}

int pt_GetOption(pt_Device *dev, PT_E_OPTION option, int *value)
{
	// trap dev == NULL or value == NULL
	if ((dev == NULL) || (value == NULL)) {
		return PT_ERR_BAD_PARAMETER;
	}

	// get option value
	switch(option) {
		case PT_OPTION_MIRROR:			// Mirror
			*value = dev->mirror;
			return PT_ERR_SUCCESS;

		case PT_OPTION_AUTOCUT:		// Auto-cutter enable/disable
			*value = dev->autocut;
			return PT_ERR_SUCCESS;

		default:
			return PT_ERR_BAD_PARAMETER;
	}
}

int pt_Print(pt_Device *dev, gdImagePtr *labels, int count)
{
	int err;
	gdImagePtr *curLabel = labels;

	// trap dev == NULL
	if (dev == NULL) {
		return PT_ERR_BAD_PARAMETER;
	}

	// trap labels == NULL
	if (labels == NULL) {
		return PT_ERR_BAD_PARAMETER;
	}

	// trap count == 0
	if (count == 0) {
		return PT_ERR_BAD_PARAMETER;
	}

	// Request current status from the printer
	if ((err = pt_GetStatus(dev)) != PT_ERR_SUCCESS) {
		return err;
	}

	// Make sure the printer has tape, and is ready
	if ((dev->errorInfo[0] != 0x00) || (dev->errorInfo[1] != 0x00)) {
		return PT_ERR_PRINTER_NOT_READY;
	}

	// Send the print initialisation commands
	//
	// ESC i M -- Set Mode
	fprintf(dev->fp, "%ciM%c", ESC,
			(dev->autocut ? 0x40 : 0x00) | (dev->mirror ? 0x80 : 0x00));

	// ESC i K -- Set Expanded Mode
	fprintf(dev->fp, "%ciK%c", ESC, 0x00);

	// ESC i R {n1} -- Set Raster Graphics Transfer Mode
	//   {n1} = 0x01 ==> Raster Graphics mode
	fprintf(dev->fp, "%ciR%c", ESC, 0x01);

	// M {n1} -- Set Compression Mode
	//   {n1} = 0x00 ==> no compression
	//   	Doesn't seem to work on the PT-2450DX...
	//   {n1} = 0x01 ==> reserved
	//   {n1} = 0x02 ==> TIFF/Packbits
	//   	But this works fine on the PT-2450DX...
	fprintf(dev->fp, "M%c", 0x02);

	// Loop over the images that were passed in
	for (int imnum=0; imnum < count; imnum++) {
		// Make sure image is the right size for this label tape
		if (gdImageSY(*curLabel) != dev->pixelWidth) {
			return PT_ERR_LABEL_TOO_WIDE;
		}

		// Trap a label with a width of zero
		// I'm not sure if this can happen, but it's a single if statement, so
		// probably worth checking anyway...
		if (gdImageSX(*curLabel) == 0) {
			return PT_ERR_LABEL_ZERO_LENGTH;
		}

		// Get the index of the colour "white" (RGB:255,255,255) in the Gd image
		int col_white = gdImageColorResolve(*curLabel, 255, 255, 255);

		// Iterate left-to-right over the source image
		for (int xpos = 0; xpos < gdImageSX(*curLabel); xpos++) {
			char bitbuf[128/8];		// 128-dot printhead, 8 bits per byte
			int rowclear = true;	// TRUE if entire row is clear, else FALSE

			// Fill bit buffer with zeroes
			memset(bitbuf, 0x00, sizeof(bitbuf));

			// Calculate left-side margin for this label size
			// Again, 128-dot printhead.
			int margin = (128 / 2) - (dev->pixelWidth / 2);

			// Copy data from the image to the bit-buffer
			for (int ypos = 0; ypos < gdImageSY(*curLabel); ypos++) {
				// Get pixel from gd, is it white?
				if (gdImageGetPixel(*curLabel, xpos, ypos) != col_white) {
					// No. Set the bit.
					int bit = 1 << (7 - ((margin+ypos) % 8));
					bitbuf[(margin+ypos) / 8] |= bit;

					// Clear the "row is clear" flag
					rowclear = false;
				}
			}

			// We now have the image data in bitbuf, and a flag to say whether
			// there are any black pixels in the buffer. We can pack full-rows
			// by sending a "Z" character, i.e. "this row is entirely blank".
			// If not, we have to send a (128/8)=16 byte chunk of image data.
			if (rowclear) {
				// Row is clear -- send a "clear row" command
				fprintf(dev->fp, "Z");
			} else {
				// Row is not clear -- send the pixel data
				//
				// TODO: the printer supports Packbits compression. Implement!
				// TODO: After Packbits is implemented, ((128/8)+1) must be
				//   changed to the length of the Packbits compressed data.
				//   (note: 128/8 is the printhead size in bytes, the +1 is for
				//   the control byte we add below...)
				fprintf(dev->fp, "G%c%c", ((128/8)+1) & 0xff, ((128/8)+1) >> 8);

				// This printer asks for Packbits compressed data. In this
				// case, we send a "run of N" control byte and fake it...
				fputc(sizeof(bitbuf) - 1, dev->fp);
				for (int i=0; i<sizeof(bitbuf); i++) {
					fputc(bitbuf[i], dev->fp);
				}
			}
		}

		// Is this the last label?
		if (imnum == (count-1)) {
			// Yes, send an End Job command
			fprintf(dev->fp, "%c", 0x1A);
		} else {
			// No, more labels to print. Send a formfeed.
			fprintf(dev->fp, "%c", 0x0C);
		}

		// Move onto the next label
		curLabel++;
	}

	// Operation successful.
	return PT_ERR_SUCCESS;
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

