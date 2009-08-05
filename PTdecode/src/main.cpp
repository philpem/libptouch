/****************************************************************************
 * ptdecode: P-touch PT-2450DX output decoder
 ****************************************************************************/

#include <cstdio>
#include <exception>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

// maximum size of a Ptouch printer head in dots
const unsigned int PT_HEAD_WIDTH = 1024;

// If defined, makes "blank row" blocks visible
//#define MAKE_BLANK_ROWS_VISIBLE

// custom exception class for file read errors
class EReadError : public exception {
	public:
		virtual const char* what() const throw()
		{
			return "Read error";
		}
};

FILE *fp;

// get next character from file
unsigned char getNext() {
	unsigned char ch;
	int i;

	i = fread(&ch, 1, 1, fp);
	if (i != 1) {
		throw EReadError();
	} else {
		return ch;
	}
}

// Handler for graphics transfer mode 1
void runGraphicsXferMode1()
{
	bool exit = false;
	unsigned int cm = -1;
	unsigned long xpos = 0;
	unsigned long ypos = 0;
	unsigned long ydim = 128;
	CImg<unsigned char> img(0, 0, 0, 0, (unsigned char)0);

	while (!exit) {
		unsigned char ch = getNext();
		unsigned int len = 0;
		unsigned int rowpos = 0;
		unsigned char row[PT_HEAD_WIDTH / 8];	// stores uncompressed row data

		switch (ch) {
			case 'M':			// Set compression mode
				ch = getNext();
				cm = ch;
				printf("Set compression mode: 0x%02X", ch);
				switch (cm) {
					case 0x02:
						printf(" (TIFF/Packbits)\n");
						break;
					default:
						printf(" *** Unknown, assuming uncompressed ***\n");
						cm = 1;
						break;
				}
				break;

			case 'Z':			// Blank raster line
				// Increment x-position and resize the image
				img.resize(xpos+1, ydim, 1, 1, 0, 0);

				// Blank the new row
				if (img.dimy() > 0) {
//					printf("Clear row: x=%lu\n", xpos);
					for (int i=0; i<img.dimy(); i++) {
#ifdef MAKE_BLANK_ROWS_VISIBLE
						img(xpos, i) = 128;
#else
						img(xpos, i) = 255;
#endif
					}
				}

				xpos++;
				break;

			case 'G':			// Graphics data row
				// decode the length
				ch = getNext();
				len = (((int)getNext()) << 8) + ch;

				// Is gfx payload compressed or uncompressed?
				if (cm == 1) {
					// Uncompressed. Read straight into the row buffer.
					while (len > 0) {
						row[rowpos++] = getNext(); len--;
					}
				} else {
					// Decompress the gfx data
					rowpos = 0;
					while (len > 0) {
						// get the prefix byte
						ch = getNext(); len--;

						// Is this a "run" (a single byte replicated) or a "copy"?
						int runlen;
						if (ch & 0x80) {
							// MSB set, it's a run
							runlen = 257 - ((int)ch);

							// Get the byte to replicate, and replicate it into the o/p buffer
							ch = getNext(); len--;
							while (runlen-- > 0) {
								row[rowpos++] = ch;
							}
						} else {
							// MSB clear, it's a copy
							runlen = ((int)ch) + 1;

							// Copy N bytes from the input stream to the output
							while (runlen-- > 0) {
								row[rowpos++] = getNext();
								len--;
							}
						}
					}
				}

				// Row decode complete. row contains the image data, and rowpos
				// contains its length in bytes. Now shuffle it into CImg...

				// If image height is less than size of image row, then make the
				// image taller.
				if (((unsigned int)img.dimy()) < (rowpos * 8)) {
					ydim = rowpos * 8;
				} else {
					ydim = img.dimy();
				}

				// Perform the Y resize if necessary, but also make Xdim=Xdim+1
				img.resize(xpos+1, ydim, 1, 1, 0, 0);

				img(xpos, ydim/2) = 128;

				// Now copy the image data...
				ypos = 0;
				for (unsigned int byte=0; byte<rowpos; byte++) {
					for (unsigned int bit=0; bit<8; bit++) {
						if (row[byte] & (0x80>>bit)) {
							img(xpos, ypos) = 0;
						} else {
							img(xpos, ypos) = 255;
						}

						// Increment y-position
						ypos++;
					}
				}

				// An entire row has been decoded. Increment x-position.
				xpos++;
				break;

			case 0x0c:		// FF
				printf("Formfeed: Print without label feed (job completed, more labels follow)\n");
				exit = true;
				break;

			case 0x1a:		// Ctrl-Z
				printf("Ctrl-Z:   Print with label feed    (job completed, no further labels)\n");
				exit = true;
				break;

			default:			// Something else
				printf("** Unrecognised command prefix in gfx mode: 0x%02x\n", ch);
				break;
		}
	}

	// Display the contents of the image
	img.display();
}

// Parse an ESC i command
void parse_esc_i()
{
	unsigned char ch = getNext();
	unsigned int tmpI;

	switch (ch) {
		case 'B':				// ESC i B: Specify baud rate
			tmpI = getNext();
			ch = getNext();
			tmpI += ((int)ch)*256;
			printf("Set baud rate:\t%d00", tmpI);
			if ((tmpI != 96) && (tmpI != 576) && (tmpI != 1152)) {
				printf(" [ILLEGAL SETTING]\n");
			} else {
				printf("\n");
			}
			break;

		case 'S':				// ESC i S: Status request
			printf("Printer status request\n");
			break;

		case 'M':				// ESC i M: Set mode
			ch = getNext();
			printf("Set mode 0x%02X:\tAutoCut %s, Mirror %s\n", ch,
					(ch & 0x40) ? "on" : "off",
					(ch & 0x80) ? "on" : "off");
			break;

		case 'd':				// ESC i d: Set margin amount (feed amount)
			tmpI = getNext();
			ch = getNext();
			tmpI += ((int)ch)*256;
			printf("Set margin:\t%d dots", tmpI);
			break;

		case 'K':				// ESC i K: Set expanded mode
			ch = getNext();
			printf("Set expanded mode 0x%02X:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n", ch,
					(ch & 0x04) ? "Half-cut on" : "Half-cut off",
					(ch & 0x08) ? "Chain-print off: last label will be fed and cut" : "Chain-print on: last label will NOT be fed or cut",
					(ch & 0x20) ? "Label end cut: when printing multiple copies, end of last label is cut" : "Label end cut off",
					(ch & 0x40) ? "High-resolution (360x720dpi)" : "Normal resolution (360x360dpi)",
					(ch & 0x80) ? "Copy-printing on (expansion buffer not cleared on form-feed)" : "Copy-printing off"
				  );
			break;

		case 'R':				// ESC i R: Set graphics transfer mode
			ch = getNext();
			printf("Set graphics transfer mode 0x%02X: ", ch);
			if (ch == 1) {
				printf("Raster graphics mode\n");
				runGraphicsXferMode1();
			} else {
				printf("\n\tUnrecognised graphics transfer mode: remainder of data may be garbage.\n");
			}
			break;

		default:
			printf("Unrecognised cmnd: ESC i 0x%02X\n", ch);
			break;
	}
}

// Parse an ESC command
void parse_esc()
{
	unsigned char ch = getNext();

	switch(ch) {
		case 'i':		// ESC i: Brother-specific extensions
			parse_esc_i();
			break;

		case '@':		// ESC @: Initialize
			printf("Initialize: clear buffer and reset print origin\n");
			break;

		default:
			printf("Unrecognised cmnd: ESC 0x%02X\n", ch);
			break;
	}
}

int main(int argc, char **argv)
{
	// check params
	if (argc != 2) {
		// wrong!
		printf("Usage: %s filename\n", argv[0]);
		return -1;
	}

	// open binary dump file
	fp = fopen(argv[1], "rb");
	if (!fp) {
		printf("Error opening source file\n");
		return -1;
	}

	try {
		while (true) {
			unsigned char ch;

			ch = getNext();

			switch (ch) {
				case 0x00:		// NULL
					printf("Null\n");
					break;
				case 0x0c:		// FF
					printf("Formfeed: Print without feed\n");
					break;
				case 0x1a:		// Ctrl-Z
					printf("Ctrl-Z:   Print with label feed\n");
					break;
				case 0x1b:		// ESC
					parse_esc();
					break;
				default:
					printf("Unrecognised cmnd: 0x%02X\n", ch);
					break;
			}
		}
	} catch (EReadError &e) {
		if (feof(fp)) {
			printf("EOF reached.\n");
		} else {
			printf("Uncaught EReadException: %s\n", e.what());
		}
	} catch (exception &e) {
		printf("Uncaught exception: %s\n", e.what());
	} catch (...) {
		printf("Uncaught and unrecognised exception. Something went *really* wrong here...\n");
	}

	// close the file
	fclose(fp);

	return 0;
}
