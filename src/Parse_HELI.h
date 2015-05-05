#ifndef __Parse_HELI__
#define __Parse_HELI__

#include "Constants.h"
#include "Globals.h"


//-----------------------------------------------------------------------------

extern TComplex heli_H1[100][EBINS];
extern TComplex heli_H2[100][EBINS];
extern TComplex heli_H3[100][EBINS];
extern TComplex heli_H4[100][EBINS];
extern Double_t heli_z[100][EBINS];
extern Double_t heli_en[EBINS];
extern Int_t heli_bin;
extern Int_t heli_nz[EBINS];

//-----------------------------------------------------------------------------

void Parse_HELI(Char_t*);
Int_t GetEnergyBin_heli();


//-----------------------------------------------------------------------------

#endif

