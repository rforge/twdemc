// to be included as part of .c file

// assumed definitions of parameter names MODEL_ID,
// PDBL_TABLE, PINT_TABLE, and PCHR_TABLE of type X(foo) X(bar)
// declares global variables
// which then can be accessed APPMID(PD)[foo]

// assumes N_XCOLS_C 2, XCOLS_CI {s,a}, N_XCOLS_N 1, XCOLS_NI {n15}

// global variables to hold parameters and the used decomposition function,
// initialized in init_soilmod_fs

// append defined MODEL_ID, 2 indirection needed to replace makro argument
#define APPMID__(a,b) a##_##b
#define APPMID_(a,b) APPMID__(a,b)
#define APPMID(a) APPMID_(a,MODEL_ID)

// possiblity to define value as a funtion of anohter macro
#define SUBST__(a) a
#define SUBST_(a) SUBST__(a)
#define SUBST(a) SUBST_(a)

// name of the initializer function
#define INIT_MOD APPMID(init_soilmod)

// integer constants to access components in global arrays
#define X(a) a,
enum { PDBL_TABLE };
enum { PINT_TABLE };
enum { PCHR_TABLE };
#undef X

// After parameters constants have been defined, we can define the global parameter variables
// Then we can access parmeter parameter kF of type double by PD[kF]
#define PD   APPMID(parms_dbl)
#define PINT APPMID(parms_int)
#define PCHR APPMID(parms_chr)

double PD[N_PDBL];		// double parameter values
int PINT[N_PINT];		// integer parameter values
cststr PCHR[N_PCHR];	// character parameter values (take care on scope where those are valid and not deallocated again)

// array of parameter names (# puts quotes around a)
#define PDBL_NAMES APPMID(parms_dbl_names)
#define PINT_NAMES APPMID(parms_int_names)
#define PCHR_NAMES APPMID(parms_chr_names)

#define X(a) #a,
const char *PDBL_NAMES[] = { PDBL_TABLE };
const char *PINT_NAMES[] = { PINT_TABLE };
const char *PCHR_NAMES[] = { PCHR_TABLE };
#undef X

// Columns for carbon and nitrogen isotopes
#define XCOLS_C APPMID(xcols_c)
#define XCOLS_N APPMID(xcols_n)
const unsigned char XCOLS_C[N_XCOLS_C] = XCOLS_CI;
const unsigned char XCOLS_N[N_XCOLS_N] = XCOLS_NI;

// units for isotopes to bring columns (XCOLS) to the same magnitude in order to decrease numerical errors
// are adjusted in init
#define IRUNITS APPMID(iRUnits)
double IRUNITS[N_XCOLS];

// define the global pointer to the decomposition function
//#define FDECS APPMID(f_decomp_S)
//double (*FDECS) (double *S, double *enz );		// pointer to decomposition function

