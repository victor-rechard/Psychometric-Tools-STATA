version 18
mata: 

real scalar lnv(real scalar theta,real scalar group,real scalar beta, real matrix thres,real rowvector resp ){
	
	real matrix DM
	
	real scalar centralterm,lastterm,res
	
	
	DM=J(rows(thres),1,0)
	for(i=1; i<=cols(thres); i++){
		DM=DM,DM[.,i]+thres[.,i]
	}
	
	lastterm=0
	centralterm=0
	for(j=1;j<=rows(thres);j++){
		centralterm=centralterm+DM[j,resp[j]+1]
		denomin=0
		for(i=0;i<=cols(thres);i++){
		denomin=denomin+exp(i*theta+(group*beta)-DM[j,i+1])
	}
	lastterm=lastterm+ln(denomin)
	}
	
	
	
	res=0.5*ln(testinfo(theta,group,beta,thres))+theta*sum(resp)-centralterm-lastterm

	return(res)
	
}

end
