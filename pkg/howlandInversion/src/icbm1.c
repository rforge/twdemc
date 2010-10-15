/* soilmod_fs.c
system("R CMD SHLIB src/soilmod_fs.c src/rc_helpers.c src/mic_c_part.c")
*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <string.h>   //memcpy
#include "rc_helpers.h"

#include "c14constants.h"


// Make sure that all global object names are unique
// Append MODEL_ID do names
// see macro APPMID and definitions in par_define.c
#define MODEL_ID icbm1

// define the constants to access rows (stateVars) and columns (isotopic components) in X
enum XROWS {Y,O,R,     N_XROWS};
enum XCOLS {c12,c14,   N_XCOLS};

// define columns in X comprising carbon or nitrogen, required in sums across elements
#define N_XCOLS_C 2
#define XCOLS_CI {c12,c14}
#define N_XCOLS_N 0
#define XCOLS_NI {}


// define constants to access components in *aux, i.e. auxiliary outputs in addition to *x
enum AUX_OUTPUT_NAMES {inputLeaf_c12,inputLeaf_c14,inputRoot_c12,inputRoot_c14,decY_c12,decY_c14,decO_c12,decO_c14,respY_c12,respY_c14,respO_c12,respO_c14,N_AUX}; //generated in ICBM1.R

// define parameters and global variables
#define PDBL_TABLE \
X(kY) X(kO) \
X(h) \
X(N_PDBL)
#define PINT_TABLE X(N_PINT)
#define PCHR_TABLE X(N_PCHR)
#include "src_incl/par_define.c"     // defines all global variables and constants based on already defined constants.

//#include "src_incl/decomp.c"		// defines the decomposition functions, assume parameter arrays and constants defined

#include "src_incl/init_soilmod.c"	// defines initializing function init_soilmod_MODEL_ID. constants  need to be declared

#define XRC(iRow,iCol) (x[(iCol)*N_XROWS+(iRow)])    //access double *x as matrix, e.g. XRC(A,n15)
#define DXRC(iRow,iCol) (dx[(iCol)*N_XROWS+(iRow)])  //access double *dx as matrix, e.g. DXRC(A,n15)

//#define SQRT_DBL_EPSILON SUBST(DBL_EPSILON)		// TODO assign sqrt(DBL_EPSILON) to global variable


// define the inputs and the initialization function. access by e.g. input[leafC12]
// note that inputs must be provided in this order as forcings argument (forcings[1] is matrix for leafC12 with columns time and value)
// note that must provide all isotopes
enum INPUT{ leafC12, rootC12, leafC14, rootC14, N_INPUT};		// matrix by column (first all components of one isotope) as in R
enum INPUT_ROWS{ leaf, root, N_INPUT_ROWS};
static double input[N_INPUT];
#define INPUTRC(iRow,iCol) (input[(iCol)*N_INPUT_ROWS+(iRow)])    //access double *input as matrix, e.g. IRC(Y,c12)
//#define import forc[0]
void APPMID(forcc) (void (* odeforcs)(int *, double *)) {
	int N=N_INPUT;  // number of variables
	odeforcs(&N, input);
}


void APPMID(deriv) (int *neq, double *t, double *x, double *dx, double *aux, int *ip) {
	if (*ip < N_AUX) error("nout should be at least N_AUX");

	//---------- check for non-finite values
	//for( unsigned short int i=0; i<(N_XROWS*N_XCOLS); i++ )
	//	if( !finite(x[i]) ) error("Derivate function provided with non-finite values.");

	//-----  loop across all element isotopes
	double decY[N_XCOLS], decO[N_XCOLS];
	for( unsigned char j=0; j<N_XCOLS; j++ ){
		decY[j] = PD[kY] * XRC(Y,j);
		decO[j] = PD[kO] * XRC(O,j);

		DXRC(Y,j) = +INPUTRC(leaf,j) +INPUTRC(root,j) -decY[j];
		DXRC(O,j) = -decO[j];		// +PD[h]*decY[j] differens between C and N, N is not respired

		aux[ inputLeaf_c12+j ] = INPUTRC(leaf,j);
		aux[ inputRoot_c12+j ] = INPUTRC(root,j);
		aux[ decY_c12+j ] = decY[j];
		aux[ decO_c12+j ] = decO[j];
	}
	DXRC(O,c14) -= XRC(O,c14)*lambda;	//radioactive decay see c14Utils

	//------ loop across carbon isotopes
	double respY[N_XCOLS_C], respO[N_XCOLS_C];
	for( unsigned char jc=0; jc<N_XCOLS_C; jc++ ){
		int j = XCOLS_C[jc];
		respY[jc] = decY[j] * (1-PD[h]);
		respO[jc] = decO[j];
		DXRC(R,j) = +respY[jc] +respO[jc];
		DXRC(O,j) += +PD[h]*decY[j];

		aux[ respY_c12+jc] = respY[jc];
		aux[ respO_c12+jc] = respO[jc];
	}

	//------ loop across nitrogen isotopes
	// O+decY


}

