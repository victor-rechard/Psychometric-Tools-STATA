program define wle ,rclass
version 18
args thres resp


mata: 

mata clear

real scalar pcmformula(real scalar theta,  real matrix thres,  real scalar resp, real scalar j){
	
	real scalar denomin
	real matrix DM

	

	DM=J(rows(thres),1,0)
	for(i=1; i<=cols(thres); i++){
		DM=DM,DM[.,i]+thres[.,i]
	}
	denomin=0
	for(i=0;i<=cols(thres);i++){
		denomin=denomin+exp(i*theta-DM[j,i+1])
	}
	return(exp(resp*theta-DM[j,resp+1])/denomin)
}



real scalar testinfo(real scalar theta , real matrix thres){
	
	real scalar res,part1,part2,proba
	
   res=0
	
	for(j=1;j<=rows(thres);j++){
			part1=0
	
	for(x=0;x<=cols(thres);x++){
		proba=pcmformula(theta,thres,x,j)
		part1=part1+proba*x^2
	}
		
		part2=0
		for(x=0;x<=cols(thres);x++){
		proba=pcmformula(theta,thres,x,j)
		part2=part2+proba*x
	}
	
	res=res+(part1-part2^2)
	
	
	}
	return(res)
}
real scalar lnv(real scalar theta, real matrix thres,real rowvector resp ){
	
	real matrix DM
	
	real scalar centralterm,lastterm,info,res
	
	
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
		denomin=denomin+exp(i*theta-DM[j,i+1])
	}
	lastterm=lastterm+ln(denomin)
	}
	
	
	
	res=0.5*ln(testinfo(theta,thres))+theta*sum(resp)-centralterm-lastterm

	return(res)
	
}



void wle_ind( todo,p,thres,resp, y, g, H){
	
	x=p
	y=lnv(x,thres,resp)
	
}

real  matrix wle(string name_thres, string names_item_res ){
	
	real rowvector res
	real matrix seuils
	real rowvector R
	transmorphic S
	seuils=st_matrix(name_thres)
	
	res=J(st_nobs(),1,.)
	
	for(i=1;i<=st_nobs();i++){
	R=st_data(i,(names_item_res))
	S=optimize_init()
	optimize_init_evaluator(S,&wle_ind())
	optimize_init_evaluatortype(S,"gf0")
	optimize_init_params(S,0)
	optimize_init_argument(S,1,seuils)
	optimize_init_argument(S,2,R)
	optimize_init_tracelevel(S,"none")
	res[i]=optimize(S)
}
return(res)

}

end

qui mata: prediction=wle("`thres'","`resp'")

qui getmata prediction
end


