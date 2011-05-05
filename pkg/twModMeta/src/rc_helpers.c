#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "rc_helpers.h"


SEXP getListElement(SEXP list, cststr str) {
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < length(list); i++)
	 if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	   elmt = VECTOR_ELT(list, i);
	   break;
	 }
	return elmt;
}


void extractParameters( const SEXP parms
	, const int nPdbl, const cststr *pDblNames, double *pd
	, const int nPint, const cststr *pIntNames, int *pi
	, const int nPchr, const cststr *pChrNames, cststr *pc
){
	// extract parameters from R list parms into arrays pd, pi and pc
	// note that referenced strings in pc may be deallocated when SEXP goes out of scope.
	// So if you need strings outside, copy the string yourself and also take care of memory deallocation.
	SEXP elmt;
	for( int i=0; i<nPdbl; i++ ){
		elmt = getListElement( parms, pDblNames[i] );
		if( elmt==R_NilValue ){
			char message[256];
			sprintf(message,"Double parameter '%s', has not been given.",pDblNames[i]);
			error(message);
		}
		pd[i] = *REAL(coerceVector(elmt, REALSXP));
	}
	for( int i=0; i<nPint; i++ ){
		elmt = getListElement( parms, pIntNames[i] );
		if( elmt==R_NilValue ){
			char message[256];
			sprintf(message,"Integer parameter '%s', has not been given.",pIntNames[i]);
			error(message);
		}
		pi[i] = *INTEGER(coerceVector(elmt, INTSXP));
	}
	for( int i=0; i<nPchr; i++ ){
		elmt = getListElement( parms, pChrNames[i] );
		if( elmt==R_NilValue ){
			char message[256];
			sprintf(message,"Character parameter '%s', has not been given.",pChrNames[i]);
			error(message);
		}
		pc[i] = (const char *) CHAR(STRING_ELT(coerceVector(elmt, STRSXP),0));
	}
}

