orderLogDen <- function(
    ### rank the log-Densities, starting with highest densities
    ss          ##<< numeric matrix ( nRec x nDen+nPar ) 
    ,nDen=1     ##<< number of densities
){
    oR <- if( nDen==1){
        order( ss[,1] )  
    }else{
        r <- apply( ss[,1:nDen], 2, rank )
        maxR <- apply(r,1,max)
        order(maxR)
    }
    ### numeric vector of indices of records with lowest rank of log-Density starting with the maximum density.
    ### For each record the maximum of the density ranks is calculated
    oR
}