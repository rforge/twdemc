typedef const char* cststr;

SEXP getListElement(SEXP list, cststr str);

void extractParameters( const SEXP parms
	, const int nPdbl, const cststr *pDblNames, double *pd
	, const int nPint, const cststr *pIntNames, int *pi
	, const int nPchr, const cststr *pChrNames, cststr *pc
);

