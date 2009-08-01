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

typedef struct {
	FILE	*fp;
	int		headMark, size, errorInfo[2];
	int		mediaWidth, mediaType, mediaLength;
	int		statusType, phaseType, phaseHi, phaseLo;
	int		notification;
} pt_Device;

// printer functions
pt_Device *pt_Initialise(char *path);
int pt_GetStatus(pt_Device *dev);
void pt_Close(pt_Device *dev);

#endif // PTOUCH_H
