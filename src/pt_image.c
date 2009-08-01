#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pt_image.h"

/**
 * Create a new pt_Image object.
 *
 * Creates a new pt_Image of a specified width and height, and initialises
 * all the pixels in its buffer to zero.
 *
 * @param	width		The width of the image.
 * @param	height		The height of the image.
 * @return	The new pt_Image, or NULL on failure.
 */
pt_Image *ptimage_Create(unsigned long width, unsigned long height)
{
	pt_Image *image;

#ifdef DEBUG
	printf("%s[%d]: image size = %u\n", __FILE__, __LINE__, sizeof(pt_Image));
	printf("%s[%d]: data  size = %u\n", __FILE__, __LINE__, sizeof(image->data[0]));
	printf("%s[%d]: create %lu x %lu image\n", __FILE__, __LINE__, width, height);
#endif

	// Allocate memory for the new image and its data buffer
	image = malloc(sizeof(pt_Image));
	if (image == NULL) return NULL;

	image->data = malloc(width * height * sizeof(image->data[0]));
	if (image->data == NULL) return NULL;

	// Set the image's parameters
	image->width = width;
	image->height = height;

	// Clear the image's data buffer to "blank"
	memset(image->data, 0, width * height * sizeof(image->data[0]));

	return image;
}

/**
 * Free a pt_Image object.
 *
 * @param	image		The image object to free.
 */
void ptimage_Free(pt_Image *image)
{
	// Make sure image is non-null
	if (image == NULL) {
		return;
	}

	// Free the image data
	if (image->data != NULL) {
		free(image->data);
	}

	// Free the image
	free(image);
}

/**
 * Get the value of a pixel in a pt_Image object.
 *
 * @param	image		Image object.
 * @param	x			X position of the pixel, zero-based.
 * @param	y			Y position of the pixel, zero-based.
 * @return	Value of the pixel, or negative on error.
 */
int ptimage_GetPixel(pt_Image *image, unsigned long x, unsigned long y)
{
	// Make sure the image is not null and that the data buffer has been
	// allocated
	if (image == NULL) {
		return -1;	// TODO: make constant
	}

	if (image->data == NULL) {
		return -1;	// TODO: make constant
	}

	// Range-check
	if ((x < 0) || (x > image->width)) {
		return -2;	// TODO: make constant
	}
	if ((y < 0) || (y > image->height)) {
		return -3;	// TODO: make constant
	}

	// Return the pixel value
	return image->data[(y*image->width)+x];
}

/**
 * Set the value of a pixel in a pt_Image object.
 *
 * @param	image		Image object.
 * @param	x			X position of the pixel, zero-based.
 * @param	y			Y position of the pixel, zero-based.
 * @param	val			New value of the pixel.
 * @return	Zero on success, or negative on error.
 */
int ptimage_SetPixel(pt_Image *image, unsigned long x, unsigned long y, unsigned char val)
{
	// Make sure the image is not null and that the data buffer has been
	// allocated
	if (image == NULL) {
		return -1;	// TODO: make constant
	}

	if (image->data == NULL) {
		return -1;	// TODO: make constant
	}

	// Range-check
	if ((x < 0) || (x > image->width)) {
		return -2;	// TODO: make constant
	}
	if ((y < 0) || (y > image->height)) {
		return -3;	// TODO: make constant
	}

	// Set the pixel value
	image->data[(y*image->width)+x] = val;

	return 0;
}

