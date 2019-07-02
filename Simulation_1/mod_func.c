#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _AdExp_reg();
extern void _izhi2007b_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," AdExp.mod");
fprintf(stderr," izhi2007b.mod");
fprintf(stderr, "\n");
    }
_AdExp_reg();
_izhi2007b_reg();
}
