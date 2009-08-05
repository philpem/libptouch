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

/// Label has a length of zero
	PT_ERR_LABEL_ZERO_LENGTH	= -4,

/// Printer is not ready
	PT_ERR_PRINTER_NOT_READY	= -5
};

/*
 * Job options
 */
typedef enum {
/// Mirror -- mirror the printed label along the long edge
	PT_OPTION_MIRROR,
/// Auto-cutter -- enable or disable automatic label cutting
	PT_OPTION_AUTOCUT
} PT_E_OPTION;

/**
 * @brief	Initialise the printer.
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
 * @brief	Close a printer device.
 *
 * Closes the connection to the printer, and destroys the pt_Device struct.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 */
void pt_Close(pt_Device *dev);

/**
 * @brief	Get the current status of the printer.
 *
 * Queries the printer for its current status, then returns the result.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 * @return	Any valid PT_ERR_* constant. The pt_Device struct passed in
 * 	is also updated with the current status of the printer.
 */
int pt_GetStatus(pt_Device *dev);

/**
 * @brief	Set a job option for the next print job.
 *
 * Sets a job option (specified by <b>option</b>) for the next print job.
 * These options include printer features like auto-cutting and mirroring
 * of the label image.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 * @param	option	One of the PT_OPTION_* constants specifying the parameter
 * 	that is to be set.
 * @param	value	The value to assign to the job option.
 * @return	<b>PT_ERR_BAD_PARAMETER:</b> Either <b>dev</b> was equal to NULL,
 * 	or the option value specified in <b>option</b> was invalid.<br>
 * 	<b>PT_ERR_SUCCESS:</b> Operation completed successfully, the current value
 * 	of the option parameter is now set to the contents of <b>value</b>.
 */
int pt_SetOption(pt_Device *dev, PT_E_OPTION option, int value);

/**
 * @brief	Get the current value of a job option for the next print job.
 *
 * Returns the current value of a job option (specified by <b>option</b>)
 * for the next print job.
 * These options include printer features like auto-cutting and mirroring
 * of the label image.
 *
 * @param	dev		A pt_Device struct created by pt_Initialise.
 * @param	option	One of the PT_OPTION_* constants specifying the parameter
 * 	that is to be returned.
 * @param	value	A pointer to an <b>int</b> that will contain the value
 * 	of the job option.
 * @return	<b>PT_ERR_BAD_PARAMETER:</b> Either <b>dev</b> or <b>value</b> was
 * 	equal to NULL, or the option value specified in <b>option</b> was invalid.<br>
 * 	<b>PT_ERR_SUCCESS:</b> Operation completed successfully, the current value
 * 	of the option parameter is now stored in <b>value</b>.
 */
int pt_GetOption(pt_Device *dev, PT_E_OPTION option, int *value);

/**
 * @brief	Print one or more labels.
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
