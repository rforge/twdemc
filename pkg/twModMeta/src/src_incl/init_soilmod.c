// assumes constants and global variables are defined
// APPMID() appends MODEL_ID
// APPMID(INIT_MOD): name of the initializing function
// (N_PDBL, APPMID(PDBL_NAMES), PD : number of parameters, global cststr array of parameter names, global array of double parameters
// (N_PDBL, APPMID(PINT_NAMES), APPMID(PINT)
// N_PCHR, APPMID(PCHR_NAMES), PCHR
// used as include after those definitions in .c file
// intializes global variables PD, pRUnits, and function pointer FDECS from examining parameter parms to R function lsoda
// cannot reuse it for different models, because then it is defined several times with different meanings.
void INIT_MOD (void(* odeparms)(int *, double *)) {
	DL_FUNC get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
	SEXP parms = get_deSolve_gparms();

	extractParameters(parms
		,N_PDBL, PDBL_NAMES, PD
		,N_PINT, PINT_NAMES, PINT
		,N_PCHR, PCHR_NAMES, PCHR
	);	// needs to check correct structure

	SEXP elmt = getListElement( parms, "modMeta" );
	if( elmt == R_NilValue) error("Parms must provide element modMeta.");
	SEXP elmIRUnits = getListElement( elmt, "iRUnits" );
	if( elmIRUnits == R_NilValue) error("modMeta must contain numeric vector iRUnits.");
	if( !isReal(elmIRUnits) ) error("iRUnits must be of type double.");
	if( LENGTH(elmIRUnits) != N_XCOLS) error("iRUnits must be of length of columns of x.");
	double *pIRUnits = REAL(coerceVector(elmIRUnits, REALSXP));
	// efficient copy to IRUNITS
	memcpy( IRUNITS, pIRUnits, N_XCOLS*sizeof(IRUNITS[0]) );

	/*
	if( strcmp(PCHR[fSDec], "schimel") == 0 ) FDECS = APPMID(decomp_schimel);
	else if( strcmp(PCHR[fSDec], "firstOrder") == 0 ) FDECS = APPMID(decomp_firstOrder);
	else if( strcmp(PCHR[fSDec], "monod") == 0 ) FDECS = APPMID(decomp_monod);
	else if( strcmp(PCHR[fSDec], "fontaine") == 0 ) FDECS = APPMID(decomp_fontaine);
	else if( strcmp(PCHR[fSDec], "product") == 0 ) FDECS = APPMID(decomp_product);
	else if( strcmp(PCHR[fSDec], "wutzler") == 0 ) FDECS = APPMID(decomp_wutzler);
	else if( strcmp(PCHR[fSDec], "min") == 0 ) FDECS = APPMID(decomp_min);
	else error("unknown identifier of decomposition function.");
	*/
}

