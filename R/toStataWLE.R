
print(seuils)
library(PP)
# Fonction toStataWLE permet d'estimer le trait latent avec WLE que stata ne peut pas faire pour l'instant:
toStataWLE_PCM = function(resp,deltahat){
  #INPUT:
  #resp: matrice des réponses aux items de tout les individus
  #deltahat: matrice des estimations des paramètres de difficultés de tout les items
  #OUTPUT:  WLE pour chaque individus 
  ##NECESSITE LE CHARGEMENT DU PACKAGE PP##
  return(PP_gpcm(resp,t(deltahat),rep(1,dim(resp)[2]),type="wle")$resPP)
}
nbitem=sum(1*grepl("item",names(data)))
resp=as.matrix(data[,paste0("item",seq(nbitem))])

data[,"Theta"]=toStataWLE_PCM(resp,seuils)