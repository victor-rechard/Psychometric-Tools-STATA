*! version 1.1 25 June 2025
*! Victor Rechard
***************
* INPUT : 
*	thres : Threshold parameters estimates as Matrix (J,M_J-1)
*	resp  : String indicates response variables
*	time  : variable time name
*	U     : Random effect estimation for each individual
*	group : variable group name
*	beta  : column vector 1 then estimations of latent regression coefficients 
*
*OUTPUT :
*	New variable named prediction with weighted likelihood estimates 
***************


program define wle_long ,rclass
version 18
args thres resp time group U beta


 mata: prediction=wle("`thres'","`resp'","`time'","`group'","`U'","`beta'")
qui getmata prediction
end


mata: 

mata clear

real scalar pcmformula(real scalar theta,real scalar ajust,real matrix thres,  real scalar resp, real scalar j){
	
	real scalar denomin
	real matrix DM

	

	DM=J(rows(thres),1,0)
	for(i=1; i<=cols(thres); i++){
		DM=DM,DM[.,i]+thres[.,i]
	}
	denomin=0
	for(i=0;i<=cols(thres);i++){
		denomin=denomin+exp(i*theta+ajust-DM[j,i+1])
	}
	return(exp(resp*theta+ajust-DM[j,resp+1])/denomin)
}



real scalar testinfo(real scalar theta ,real scalar ajust, real matrix thres){
	
	real scalar res,part1,part2,proba
	
   res=0
	
	for(j=1;j<=rows(thres);j++){
			part1=0
	
	for(x=0;x<=cols(thres);x++){
		proba=pcmformula(theta,ajust,thres,x,j)
		part1=part1+proba*x^2
	}
		
		part2=0
		for(x=0;x<=cols(thres);x++){
		proba=pcmformula(theta,ajust,thres,x,j)
		part2=part2+proba*x
	}
	
	res=res+(part1-part2^2)
	
	
	}
	return(res)
}
real scalar lnv(real scalar theta,real scalar ajust, real matrix thres,real rowvector resp ){
	
	real matrix DM
	
	real scalar centralterm,lastterm,res
	
	
	DM=J(rows(thres),1,0)
	for(i=1; i<=cols(thres); i++){
		DM=DM,DM[.,i]+thres[.,i]
	}
	
	lastterm=0
	centralterm=0
	for(j=1;j<=rows(thres);j++){
		if (resp[j]!=.){
		centralterm=centralterm+DM[j,resp[j]+1]
		}
		denomin=0
		for(i=0;i<=cols(thres);i++){
		denomin=denomin+exp(i*theta+ajust-DM[j,i+1])
	}
	if (resp[j]!=.){
	lastterm=lastterm+ln(denomin)
	}
	}
	
	
	
	res=0.5*ln(testinfo(theta,ajust,thres))+theta*sum(resp)-centralterm-lastterm+rows(thres)*ajust

	return(res)
	
}






void wle_ind( todo,x,thres,resp,ajust, y, g, H){
	
	y=lnv(x,ajust,thres,resp)
	
}


real  matrix wle(string name_thres,  string names_item_res,  string name_group, string name_time, string name_U ,string name_beta ){
	
	real rowvector res ,U,G,T
	real matrix seuils,A,Ajust
	transmorphic S
	seuils=st_matrix(name_thres)
	beta=st_matrix(name_beta)
	res=J(st_nobs(),1,.)
	
	U=st_data(.,(name_U))
	G=st_data(.,(name_group))
	T=st_data(.,(name_time))
	
	A=U,G,T,G:*T

	Ajust=A*beta
	
	
	for(i=1;i<=st_nobs();i++){
	R=st_data(i,(names_item_res))
	if(missing(R)<rows(seuils)){
	S=optimize_init()
	optimize_init_evaluator(S,&wle_ind())
	optimize_init_evaluatortype(S,"gf0")
	optimize_init_params(S,0)
	optimize_init_argument(S,1,seuils)
	optimize_init_argument(S,2,R)
	optimize_init_argument(S,3,Ajust[i])
	optimize_init_tracelevel(S,"none")
	res[i]=optimize(S)
	}
	else {
		res[i]=.
	}
}
return(res)

}

end
