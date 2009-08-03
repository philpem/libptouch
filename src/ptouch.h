/****************************************************************************
 * Project:   P-touch printer driver library
 * Developer: Philip Pemberton
 * Purpose:   Make Brother P-touch (PT-series) printers do something besides
 *            gather dust.
 *
 *            Currently supports:
 *              PT-2450DX
 ****************************************************************************/

#ifndef PTOUCH_H
#define PTOUCH_H

#include <gd.h>

/**
 * @brief	Device information structure
 *
 * This is used to store the state of the printer as of the last call to
 * pt_GetStatus(), the current job settings, and other details required
 * by the printer driver.
 */
typedef struct {
	/// Reference to the printer device
	FILE	*fp;
	/// Error information
	int		errorInfo[2];
	/// Label width, type and length
	int		mediaWidth, mediaType, mediaLength;
	/// Label width in pixels
	int		pixelWidth;
	/// Status type, phase type, and phase number
	int		statusType, phaseType, phaseHi, phaseLo;
	/// Notification number
	int		notification;

	/// Print parameter: autocutter enable
	int		autocut;
	/// Print parameter: mirror printing enable
	int		mirror;
	/// Print parameter: half-cut enable
	/// Print parameter: chainprint enable
	/// Print parameter: label end cut
} pt_Device;

/*
 * Function return codes
 */
enum {
/// Operation completed successfully
	PT_ERR_SUCCESS				= 0,

/// Data transfer timed out
	PT_ERR_TIMEOUT				= -1,

/// Invalid parameter
	PT_ERR_BAD_PARAMETER		= -2,

/// Label image is too large for this label tape
	PT_ERR_LABEL_TOO_WIDE		= -3,

/// Printer is not ready
	PT_ERR_PRINTER_NOT_READY	= -4
};


/**
 * @brief	Initialise the printer
 * 
 * Initialises the printer and returns a pointer to a pt_Device struct
 * describing it.
 *
 * @param	path	Path to the printer device file (e.g. /dev/usb/lp0)
 * @return	On success, a pt_Device struct referring to the printer.
 * 			On failure, NULL.
 */
pt_Device *pt_Initialise(char *path);

/**
 * @brief	Close a printer device
 *
 * Closes the connection to the printer, and destroys the pt_Device struct.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 */
void pt_Close(pt_Device *dev);

/**
 * @brief	Get the current status of the printer
 *
 * Queries the printer for its current status, then returns the result.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 * @return	Any valid PT_ERR_* constant. The pt_Device struct passed in
 * 	is also updated with the current status of the printer.
 */
int pt_GetStatus(pt_Device *dev);

/**
 * @brief	Print one or more labels
 *
 * Takes a pointer to an array of Libgd images, and prints each of them.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 * @param	labels	A pointer to an array of gdImagePtr objects containing
 *	the labels to be printed.
 * @param	count	The number of labels to be printed.
 * @return	Any valid PT_ERR_* constant.
 */
int pt_Print(pt_Device *dev, gdImagePtr *labels, int count);

#endif // PTOUCH_H
