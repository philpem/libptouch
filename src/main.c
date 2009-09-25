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
#include <gd.h>
#include <gdfonts.h>
#include <gdfontl.h>
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

	// Get printer status
	int err;
	if ((err = pt_GetStatus(dev)) != PT_ERR_SUCCESS) {
		printf("getstatus error %d\n", err);
		return -1;
	}

	// Figure out image width from DPI
	// Multiply by 10, that way we can use integer math instead if floating
	// point for the DPI/width calculations.
	unsigned long labelLength = ((20.0 * dev->dpiLabel * 10) / 254) - 2;

	// quick test
	printf("label of 256 dots = %lf in\n", (256.0 / dev->dpiLabel));
	printf("label of 256 dots = %lf mm\n", (256.0 / dev->dpiLabel) * 25.4);
	printf("label of %ld dots = %lf mm\n", labelLength, (((double)labelLength) / dev->dpiLabel) * 25.4);

	// Create a label
	gdImagePtr im;
	int white, black;
	im = gdImageCreate(labelLength, dev->pixelWidth);
	white = gdImageColorAllocate(im, 255, 255, 255);
	black = gdImageColorAllocate(im, 0, 0, 0);
/*	gdImageString(im, gdFontGetLarge(), 0, 0, "!!123!! Test label !!ABC!!", black);
	gdImageString(im, gdFontGetLarge(), 10, 10, "!!123!! Test label !!ABC!!", black);
	gdImageString(im, gdFontGetLarge(), 20, 20, "!!123!! Test label !!ABC!!", black);
	gdImageString(im, gdFontGetLarge(), 30, 30, "!!123!! Test label !!ABC!!", black);
	gdImageString(im, gdFontGetLarge(), 40, 40, "!!123!! Test label !!ABC!!", black);
*/

//#define FONT_BLUEHIGHWAY

#if defined(FONT_BLUEHIGHWAY)
	char *font_normal = "./bluehighway/bluehigh.ttf";
	char *font_bold   = "./bluehighway/bluebold.ttf";
#else
	char *font_normal = "./Arial.ttf";
	char *font_bold   = "./Arial_Bold.ttf";
#endif

	char *pn = "SMD-001/01";
	char *strings[] = {
		"Surface mount resistor",
		"0805, \u00BCW (SMD-001/01)",
		"$33\u03A9",
	};

	int curtop = 3;
	int brect[8];
	char *fte;

	gdFTStringExtra extra;
	extra.flags = gdFTEX_RESOLUTION | gdFTEX_DISABLE_KERNING;
	extra.hdpi = dev->dpiLabel;
	extra.vdpi = dev->dpiPrinthead;

	//// PART NUMBER
	// calc bounding box
	fte = gdImageStringFTEx(NULL, &brect[0], 0, font_bold, 9.0, 0.0, 0, 0, pn, &extra);
	if (fte != NULL) printf("freetype error: %s\n", fte);

	// draw P/N
	fte = gdImageStringFTEx(im, NULL, 0-black, font_bold, 9.0, 0.0, 3+abs(brect[6]), curtop+abs(brect[7]), pn, &extra);
	if (fte != NULL) printf("freetype error: %s\n", fte);

	// update Y position
	curtop += abs(brect[7]) + 2;

	//// draw description
	for (int i=0; i<(sizeof(strings)/sizeof(strings[0])); i++) {
		char *str = strings[i];
		char *font = font_normal;
		float fontsize = 7.0;

		if (*str == '$') { font = font_bold; fontsize = 8.0; str++; }

		fte = gdImageStringFTEx(NULL, &brect[0], 0, font, fontsize, 0.0, 0, 0, str, &extra);
		if (fte != NULL) printf("freetype error: %s\n", fte);

		fte = gdImageStringFTEx(im, NULL, 0-black, font, fontsize, 0.0, 3+abs(brect[6]), curtop+abs(brect[7]), str, &extra);
		if (fte != NULL) printf("freetype error: %s\n", fte);

		curtop += abs(brect[7]) + 1;
	}

	// dump the image (for testing purposes)
	FILE *fp = fopen("labeldump.png", "wb");
	gdImagePng(im, fp);
	fclose(fp);

	// Set job options
	pt_SetOption(dev, PT_OPTION_AUTOCUT, 1);
	pt_SetOption(dev, PT_OPTION_MIRROR, 1);
	pt_SetOption(dev, PT_OPTION_SEPARATOR, 1);

	// Print the label
	printf("Print state code: %d\n", pt_Print(dev, &im, 1));

	gdImageDestroy(im);

	// Close the printer device
	pt_Close(dev);

	return 0;
}
