

time = c(3,4,5,6,7,8,9,10,11,12,13,14,15)
infected = c(31,82,216,299,269,242,190,125,81,52,25,22,7)
data=cbind(time, infected)

library(deSolve) #install.packages("deSolve")
library(truncnorm) #install.packages("truncnorm")
library(gridExtra) #install.packages("gridExtra")
library(cowplot) #install.packages("cowplot")

##SIR model
SIR<-function(t,x,parms){
  ##taille de chaque compartiment et de la population 
  S = x[1]
  I = x[2]
  R = x[3]
  Z = x[4]
  N = x[1]+x[2]+x[3]
  ##valeurs des parametres
  beta = parms["beta"]
  gamma = parms["gamma"]
  ##variations
  dS=-beta*S*I/N 
  dI=beta*S*I/N-gamma*I 
  dR=gamma*I 
  dZ=beta*S*I/N
  res = c(dS,dI,dR,dZ)
  list(res)
}


simulate_SIR=function(parameters){
  #parameters
  parms = c(parameters["beta"],parameters["gamma"])
  N=parameters["N"]
  #initial conditions
  init <- c(N-parameters["initI"],parameters["initI"],0,0)
  #simulation
  temps <- seq(0,15)
  solveSIR <- lsoda(y =init, times=temps, func = SIR,
                    parms = parms) 
  solutionSIR=as.data.frame(solveSIR)
  names(solutionSIR)=c("time","S","I","R","Z") 
  #merge with data
  sir_data=merge(data,solutionSIR) 
  return(sir_data)
}

####################################################################################################
####################################################################################################
#############################################question 1#############################################
####################################################################################################
####################################################################################################



b=1.7 # taux de transmission
g=0.44 #taux de "guérison"
test_parameter=function(b,g){
  theta_init =c("beta"=b,"gamma"=g,"initI"=1, "N"=763) 
  simul= data.frame(simulate_SIR(theta_init))
  #first graph : 
  df <- data.frame('time' = time, 'real_infected' = infected, 'simulated'=simul$I)
  plot(df$time, df$real_infected, xlim = c(0, 18), ylim = c(0, 300), type = "l", col = "red", xlab = "Time (periods)", ylab = "Number of infected", main = "Infection vs predictions")
  lines(df$time, df$simulated, col = "blue")}

test_parameter(1.7,0.44)

#calcul du R_0 ! 

R_0 <- function(b,g){
  R_0 = b/g
  return(R_0)
}

#on produit une pluralite de combinaison de b dans 0 jusqu'a 2 et g dans 0 jusquà 1. 
#et on graphe cela dans un plot 3d. 
combines <- matrix(nrow = length(seq(0,2,by=0.1))*length(seq(0,1,by=0.05)),ncol=3)
i=1
for (b in seq(0.1,2,by=0.1)){
  for (g in seq(0.05,1,by=0.05)){
    combines[i,1] = b
    combines[i,2] = g
    combines[i,3] <- R_0(b,g)
    i= i +1
  }
}


library(plotly)
fig <- plot_ly(x = combines[,1], y = combines[,2], z = combines[,3],names= c('beta','gamma','R_0'))
fig

#trois courbe possible 

#(i) R_0 = 1 <=> gamma = beta 
test_parameter(0.5,0.5)

#(ii) R_0 < 1 <=> beta < gamma : l'epidemie ne décolle pas 
test_parameter(0.7,0.9)

#(iii) R_0 > 1 <=> beta > gamma : l'epidemie explose puis descend rapidement car l'entierete de la population est infecte tres vite 
test_parameter(1.8,0.4)


####################################################################################################
####################################################################################################
#############################################question 2#############################################
####################################################################################################
####################################################################################################


#question a) 
#vraisemblance en la loi de poisson (lambda = X_sim)
#fonction de densité à priori (Unif[0:10] pour beta et Unif[0:1] pour gamma)
#question b) 
#distribution à posteri = distribution à priori * vraisemblance d'une poisson (X_sim) / cst de normalisation 
#question c) 
#Calcul de la vraisemblance 

vraisemblance <- function(b,g){
  theta_init =c("beta"=b,"gamma"=g,"initI"=1, "N"=763) 
  simul= data.frame(simulate_SIR(theta_init))
  df <- data.frame('time' = time, 'real_infected' = infected, 'simulated'=simul$I)
  value_vr <- prod(dpois(df$real_infected,df$simulated))
  return(value_vr)}
#ici avec b = 1.7 et g = 0.44 

vraisemblance(1.7,0.44) #2.636124e-25 proche de  0 

#test avec un jeu de parametre mauvais : b =g  = 0.7 

vraisemblance(0.7,0.7) # = 0 (inf à 2.636124e-25) (c'est logique mais la difference est tres faible)

#On calcule la distribution a posteriori en utilisant les fonctions dunif pour les densites a priori des parametres et et la fonction de vraisemblance au dessus ! 
#on ignore la constante de normalisation en temps normal c'est une valeur à une constante pres

distrib_post <-function(b,g){
  post <- dunif(b,min=0,max=10)*dunif(g,min=0,max=1)*vraisemblance(b,g)
  return(post)
}

#avec notre jeu de parametre : 2.636124e-26
distrib_post(1.7,0.44)

####################################################################################################
####################################################################################################
#############################################question 3#############################################
####################################################################################################
####################################################################################################

#Algorithme MCMC (Metropolis-Hasting)
#initiation b= 1.5 et g = 0.5 et N = 5000 (nb de simulation)

test_parameter(1.5,0.5) #on est pas mal en apparence (on devrait converger vite)
#on va construire une matrice de dimension (Nx3), troisieme colonne etant le postiriori retenu


#on tire dans une loi normale avec valeur absolue pour beta et normale tronque pour g 
#voir bonus => choix de sigma s1/s2 
loi_instrumentales <- function(b,g,s1,s2){
  g_c <- rtruncnorm(1, a=0, b=1, mean = g, sd = s2)
  b_c <- rnorm(1,mean = abs(b),sd = s1)
  return(c(b_c,g_c))
}


MCMC_metropolis_hasting <- function(b,g,N,s1,s2){
  #initiation
  thetas <- matrix(nrow = N+1, ncol = 4)
  thetas[1,1] <- b 
  thetas[1,2] <- g 
  thetas[1,3] <- distrib_post(b,g)
  thetas[1,4] <- 0
  #on lance la boucle ! 
  for (i in 1:N){
    val_b <- thetas[i,1]
    val_g <- thetas[i,2]
    candits <- loi_instrumentales(val_b,val_g,s1,s2)
    candidat_b <-  candits[1]
    candidat_g <- candits[2]
    post_cand <- distrib_post(candidat_b,candidat_g)
    post_avant <- thetas[i,3]
    r <- post_cand/post_avant
    if (r == 'NA'){
      print('sigma trop fort')
    }
    else{
      if (r >1){
        thetas[i+1,1] <- candidat_b
        thetas[i+1,2] <- candidat_g
        thetas[i+1,3] <- post_cand
        thetas[i+1,4] <- 1
      }
      else {
        seuil  <- runif(1)
        if (seuil > r){
          thetas[i+1,1] <- thetas[i,1]
          thetas[i+1,2] <- thetas[i,2]
          thetas[i+1,3] <- thetas[i,3]
          thetas[i+1,4] <- 0
        }
        else{
          thetas[i+1,1] <- candidat_b
          thetas[i+1,2] <- candidat_g 
          thetas[i+1,3] <- post_cand
          thetas[i+1,4] <- 1
        }
      } }}
  return(thetas)
}

thetas = MCMC_metropolis_hasting(1.5,0.5,5000,0.3,0.09)

#on peut jouer aussi bien sur les parametres de bases que sur la valeur des variances des loi instrumentales 

#show the best parameters => b = 1.66, g = 0.453
thetas[N,]
#let's plot it 
test_parameter(thetas[N,1],thetas[N,2])

#plot the result: 

#beta ? 

#chaines de markov convergentes ! 
plot(thetas[,1],type='l')

hist(thetas[,1])

#gamma 

plot(thetas[,2],type='l')

hist(thetas[,2])

####################################################################################################
####################################################################################################
############################################ Bonus : variances optimales #############################################
####################################################################################################
####################################################################################################


#on cherche les sigmas optimaux :
# pour cela on va enregistrer le taux de succes dans une nouvelle variable
# on va simuler plusieurs sigma different pour trouver la valeur ou le taux de succes et egale à 1/4 (0.25)

#calcul proba de succes card(succès) / card(echec)

proba_succes <- function(thetas){
  succes <- table(thetas[,4])[2]
  echec <- table(thetas[,4])[1]
  return(succes/echec)
}

p = proba_succes(thetas)

#on lance des combinaisons différentes de sigma (sigma_beta de 0.05 jusqu a 0.5 et sigma_gamma de 0.05 jusqu a 0.2)

variance_boucle <- function(seq_b,seq_g,N){
  sigma_opti <- matrix(nrow = length(seq_b)*length(seq_g), ncol = 3)
  i = 1 
  for (s1 in seq_b){
    for (s2 in seq_g){
      print(paste("sigma_b =",as.character(s1),'and sigma_g = ',as.character(s2)))
      sigma_opti[i,1] = s2
      sigma_opti[i,2] = s1
      thetas = MCMC_metropolis_hasting(1.5,0.5,N,s1,s2)
      p = proba_succes(thetas)
      sigma_opti[i,3] = p
      i = i +1 
    }}
  plot(sigma_opti[,3],type='l')
  abline(h=0.25,lty=2,lwd=1,color='purple')
}

#long aussi 
#variance_boucle(seq(0.05,0.5,by=0.05),seq(0.05,2,by=0.05),100)


#pique dans les petites valeurs => on augmente le nombre de simulation N = 2000 
# et on passe de sigma_g de 0.01 - 0.08 et sigma_b de 0.01 - 0.2 (très très long)

#CODE TRES LONG, les graphiques sont dans le latex 

#variance_boucle(seq(0.01,0.2,by=0.005),seq(0.01,0.2,by=0.005),2000)


#c est entre les combinaisons s1 in [0.005:0.025] et s2 in [0.05:0.08]

#val <- sigma_opti[400:550,3]>0.25 
#variance_optimales <- sigma_opti[431,] #bingo 

#sigma optimaux => s1 = 0.015 et s2 = 0.065

#on relance notre algorithme avec ces paramètres => 

thetas = MCMC_metropolis_hasting(1.5,0.5,5000,0.015,0.065)

#beta opti = 1.66
plot(thetas[,1],type='l')
#Burn in jusqu a la période 400 à peu près, on prend la moyenne comme b* 
N = 5000
plot(thetas[400:N,1],type='l')
abline(h=mean(thetas[400:N,1]),color="red")

print(paste("beta opti =" ,as.character(mean(thetas[400:N,1]))))

hist(thetas[,1])

#gamma opti = 0.448

plot(thetas[,2],type='l')
#Burn in jusqu a la période 400 à peu près aussi, on prend la moyenne comme g* 
plot(thetas[400:N,2],type='l')
abline(h=mean(thetas[400:N,2]),color="red")
print(paste("gamma opti =" ,as.character(mean(thetas[400:N,2]))))
hist(thetas[,2])

#le modèle parfait selon l'approche bayésienne ? 

test_parameter(1.66,0.448) #estimation bayésienne de theta 

test_parameter(1.5,0.5) #répartion initial de theta 


####################################################################################################
####################################################################################################
############################################ Question 4 : Modele SEIR  #############################################
####################################################################################################
####################################################################################################


#introduction d'une duree d incubation dans le modele SIR => SEIR

#Je reprend le code de l'enonce : 

SEIR<-function(t,x,parms){
  ##taille de chaque compartiment et de la population
  S = x[1]
  E = x[2]
  I = x[3]
  R = x[4]
  Z = x[5]
  N = x[1]+x[2]+x[3]+x[4]
  ##valeurs des parametres
  beta = parms["beta"]
  gamma = parms["gamma"]
  zeta = parms["zeta"]
  ##variations
  dS=-beta*S*I/N
  dE=beta*S*I/N-zeta*E
  dI=zeta*E-gamma*I
  dR=gamma*I
  dZ=beta*S*I/N
  res = c(dS,dE,dI,dR,dZ)
  list(res)
}

simulate_SEIR<-function(parameters){
  #parameters
  parms = c(parameters["beta"],parameters["gamma"],parameters["zeta"])
  N=parameters["N"]
  #initial conditions
  init <- c(N-parameters["initI"],0,parameters["initI"],0,0)
  #simulation
  temps <- seq(0,50)
  solveSEIR <- lsoda(y =init, times=temps, func = SEIR,
                     parms = parms)
  solutionSEIR=as.data.frame(solveSEIR)
  names(solutionSEIR)=c("time","S","E","I","R","Z")
  #merge with data
  seir_data=merge(data,solutionSEIR,all=TRUE)
  return(seir_data)
}
theta_init =c("beta"=1.66,"gamma"=0.448,"zeta"=0.5,"initI"=1, "N"=763)
simul=simulate_SEIR(theta_init)


#on augmente le nombre de date pour pas avoir de valeurs manquantes 
time = c(0:50)
#on introduit le nouvel etat (E) dans le nombre d infecte et le nombre de periode
#logique vu que chaque infecte passe 2 jours dans l etat E avant de rentrer dans I.
infected = c(c(31,82,216,299,269,242,190,125,81,52,25,22,7),rep(c(0),38))
data=cbind(time, infected)


test_parameter_SEIR=function(b,g,z){
  theta_init =c("beta"=b,"gamma"=g,'zeta'=z,"initI"=1, "N"=763) 
  theta_init_SIR =c("beta"=b,"gamma"=g,"initI"=1, "N"=763) 
  simul= data.frame(simulate_SEIR(theta_init))
  simul2= data.frame(simulate_SIR(theta_init_SIR))
  simul['I_S'] = NA
  simul['I_S'][1:16,] = simul2$I
  df <- data.frame('time' = time,'real_infected' = infected, 'simulated_SIER'=simul$I, "simulated_SIR" = simul$I_S)
  plot(df$time, df$real_infected, xlim = c(0, 50), ylim = c(0, 300), type = "l", col = "red", xlab = "Time (periods)", ylab = "Number of infected", main = "Infection vs predictions")
  lines(df$time, df$simulated_SIER, col = "blue")
  lines(df$time, df$simulated_SIR, col = "green")
}
test_parameter_SEIR(1.66,0.448,0.5)

#on remarque que l introduction de zeta, le pique de l'epidemie est bien plus lointain et est moins intense car il est dillue par la guerison qui intervient durant la latence. 



\end{lstlisting}
\end{document}