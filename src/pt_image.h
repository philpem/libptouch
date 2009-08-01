#ifndef PT_IMAGE_H
#define PT_IMAGE_H

/**
 * A storage class for simple images.
 *
 * Data is stored as bytes, as this makes the math easier (and modern
 * computer systems have plenty of RAM anyway).
 */
typedef struct {
	unsigned long	width, height;
	unsigned char	*data;
} pt_Image;

pt_Image *ptimage_Create(unsigned long width, unsigned long height);
void ptimage_Free(pt_Image *image);
int ptimage_GetPixel(pt_Image *image, unsigned long x, unsigned long y);
int ptimage_SetPixel(pt_Image *image, unsigned long x, unsigned long y, unsigned char val);

#endif	// PT_IMAGE_H
