#include <R_ext/RS.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .Fortran calls */
extern void F77_NAME(bdat20)(void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *H, void *Hx, void *Hkz, void *Skz, void *Az, void *Hsh, void *Zsh, void *Zab, void *Sokz, void *Skl, void *Vol, void *LDSort, void *Bhd, void *Ifeh, void *FixLngDef, void *NMaxFixLng, void *FixLng, void *NFixLng);
extern void F77_NAME(vbdat20)(void *n, void *vBDATBArtNr, void *vD1, void *vH1, void *vD2, void *vH2, void *vH, void *vHx, void *vHkz, void *vSkz, void *vAz, void *vHsh, void *vZsh, void *vZab, void *vSokz, void *vSkl, void *vVol, void *vLDSort, void *vBhd, void *vIfeh, void *vFixLngDef, void *vNMaxFixLng, void *vFixLng, void *vNFixLng);
extern void F77_NAME(bdatd2h2trans)(void *wBDATBArtNr, void *wD1, void *wH1, void *wD2, void *wH2, void *wHges);
extern void F77_NAME(bdatdmrhx)(void *wBDATBArtNr, void *wD1, void *wH1, void *wD2, void *wH2, void *wHges, void *wHx, void *wIErr, void *wDmRHx);
extern void F77_NAME(vbdatdmrhx)(void *n, void *vBDATBArtNr, void *vD1, void *vH1, void *vD2, void *vH2, void *vHges, void *vHx, void *vIErr, void *vDmRHx);
extern void F77_NAME(bdatdorhx)(void *wBDATBArtNr, void *wD1, void *wH1, void *wD2, void *wH2, void *wHges, void *wHx, void *wIErr, void *wDoRHx);
extern void F77_NAME(vbdatdorhx)(void *n, void *vBDATBArtNr, void *vD1, void *vH1, void *vD2, void *vH2, void *vHges, void *vHx, void *vIErr, void *vDoRHx);
extern void F77_NAME(bdatvoldhmr)(void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *Hges, void *DHGrz, void *HDHGrz, void *SekLng, void *IErr, void *VolDHmR);
extern void F77_NAME(bdatvoldhor)(void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *Hges, void *DHGrz, void *HDHGrz, void *SekLng, void *IErr, void *VolDHoR);
extern void F77_NAME(bdatvolabmr)(void *wBDATBArtNr, void *wD1, void *wH1, void *wD2, void *wH2, void *wHges, void *wA, void *wB, void *wSekLng, void *wIErr, void *wVolABmr);
extern void F77_NAME(vbdatvolabmr)(void *n, void *vBDATBArtNr, void *vD1, void *vH1, void *vD2, void *vH2, void *vHges, void *vA, void *vB, void *vSekLng, void *vIErr, void *vVolABmr);
extern void F77_NAME(bdatvolabor)(void *wBDATBArtNr, void *wD1, void *wH1, void *wD2, void *wH2, void *wHges, void *wA, void *wB, void *wSekLng, void *wIErr, void *wVolABor);
extern void F77_NAME(vbdatvolabor)(void *n, void *vBDATBArtNr, void *vD1, void *vH1, void *vD2, void *vH2, void *vHges, void *vA, void *vB, void *vSekLng, void *vIErr, void *vVolABor);
extern void F77_NAME(bdathxdx)(void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *H, void *Dx, void *Hx, void *IErr);
extern void F77_NAME(vbdathxdx)(void *n, void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *H, void *Dx, void *Hx, void *IErr);
extern void F77_NAME(bdathxdxor)(void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *H, void *Dx, void *Hx, void *IErr);
extern void F77_NAME(vbdathxdxor)(void *n, void *BDATBArtNr, void *D1, void *H1, void *D2, void *H2, void *H, void *Dx, void *Hx, void *IErr);
extern void F77_NAME(bdatrinde2hx)(void *wBDATBArtNr, void *wD1, void *wH1, void *wD2, void *wH2, void *wHges, void *wHx, void *wIErr, void *wRinde2Hx);
extern void F77_NAME(vbdatrinde2hx)(void *n, void *vBDATBArtNr, void *vD1, void *vH1, void *vD2, void *vH2, void *vHges, void *vHx, void *vIErr, void *vRinde2Hx);
extern void F77_NAME(biomasse)(void *BdatBart, void *D13, void *D2, void *H2, void *H, void *Biom);
extern void F77_NAME(vbiomasse)(void *n, void *vBdatBart, void *vD13, void *vD2, void *vH2, void *vH, void *vBiom);
extern void F77_NAME(bdatformtarif)(void *Tarif, void *BDATBArtNr, void *D, void *H, void *MwQ03BWI);
extern void F77_NAME(vbdatformtarif)(void *n, void *vTarif, void *vBDATBArtNr, void *vD, void *vH, void *vMwQ03BWI);

static const R_FortranMethodDef FortranEntries[] = {
  {"bdat20",          (DL_FUNC) &F77_NAME(bdat20),        23},
  {"vbdat20",         (DL_FUNC) &F77_NAME(vbdat20),       24},
  {"bdatd2h2trans",   (DL_FUNC) &F77_NAME(bdatd2h2trans),  6},
  {"bdatdmrhx",       (DL_FUNC) &F77_NAME(bdatdmrhx),      9},
  {"vbdatdmrhx",      (DL_FUNC) &F77_NAME(vbdatdmrhx),    10},
  {"bdatdorhx",       (DL_FUNC) &F77_NAME(bdatdorhx),      9},
  {"vbdatdorhx",      (DL_FUNC) &F77_NAME(vbdatdorhx),    10},
  {"bdatvoldhmr",     (DL_FUNC) &F77_NAME(bdatvoldhmr),   11},
  {"bdatvoldhor",     (DL_FUNC) &F77_NAME(bdatvoldhor),   11},
  {"bdatvolabmr",     (DL_FUNC) &F77_NAME(bdatvolabmr),   11},
  {"vbdatvolabmr",    (DL_FUNC) &F77_NAME(vbdatvolabmr),  12},
  {"bdatvolabor",     (DL_FUNC) &F77_NAME(bdatvolabor),   11},
  {"vbdatvolabor",    (DL_FUNC) &F77_NAME(vbdatvolabor),  12},
  {"bdathxdx",        (DL_FUNC) &F77_NAME(bdathxdx),       9},
  {"vbdathxdx",       (DL_FUNC) &F77_NAME(vbdathxdx),     10},
  {"bdathxdxor",      (DL_FUNC) &F77_NAME(bdathxdxor),     9},
  {"vbdathxdxor",     (DL_FUNC) &F77_NAME(vbdathxdxor),   10},
  {"bdatrinde2hx",    (DL_FUNC) &F77_NAME(bdatrinde2hx),   9},
  {"vbdatrinde2hx",   (DL_FUNC) &F77_NAME(vbdatrinde2hx), 10},
  {"biomasse",        (DL_FUNC) &F77_NAME(biomasse),       6},
  {"vbiomasse",       (DL_FUNC) &F77_NAME(vbiomasse),      7},
  {"bdatformtarif",   (DL_FUNC) &F77_NAME(bdatformtarif),  5},
  {"vbdatformtarif",  (DL_FUNC) &F77_NAME(vbdatformtarif), 6},
  {NULL, NULL, 0}
};

void R_init_rBDAT(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
