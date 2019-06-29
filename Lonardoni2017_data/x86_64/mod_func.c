#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _AdExp_reg(void);
extern void _pregen_reg(void);
extern void _stpsyn_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," nrnMod//AdExp.mod");
    fprintf(stderr," nrnMod//pregen.mod");
    fprintf(stderr," nrnMod//stpsyn.mod");
    fprintf(stderr, "\n");
  }
  _AdExp_reg();
  _pregen_reg();
  _stpsyn_reg();
}
