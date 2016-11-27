library(foreign)
library(lattice)
library(rstan)
library(plyr)
library(MASS)
library(systemfit)
library(stargazer)
library(xtable)
rstan_options(auto_write = TRUE)
options(mc.cores = 35)

rm(list=ls())

setwd("~/Dropbox/CospAnalysis/Analysis/")

######## Read old cospdata
fulldb <- read.dta("CollapsedCospDB.dta")
kanthak <- read.csv("KanthakCongressData.csv")
kanthak <-kanthak[,c(4:6,
                     13:15)]
kanthak <- cbind(kanthak,kanthak[,5]+kanthak[,6])
kanthak <- kanthak[,-c(5,6)]
colnames(kanthak) <- c("name","state","cd", "congress","leadership")
fulldb <- merge(fulldb,kanthak,all=TRUE)

fulldb$nonAlone <- fulldb$nSponsoredBills - fulldb$TotAlone
fulldb <- subset(fulldb,apply(fulldb,1,function(x){!all(c(x[c(23:25)])==0)}))

### Example plots for the beginning of the paper
fulldb$Loner <- with(fulldb,TotAlone/nSponsoredBills)
fulldb$sh.nalone <- 1-fulldb$Loner
fulldb$Partisan <- with(fulldb,(sh.nalone)*TotCopartisans/(TotCopartisans+TotOpposition))
fulldb$Dissident <- 1-(fulldb$Partisan+fulldb$Loner)
real.simplex <- with(fulldb,data.frame(Partisan=Partisan,Dissident=Dissident,Loner=Loner))
real.simplex$group <- "All observations"
young.ind <- which(fulldb$congress==94&fulldb$name=="YOUNGCWBILL")
real.simplex <- rbind(real.simplex,real.simplex[young.ind,])
real.simplex[nrow(real.simplex),"group"] <- "Single portfolio"
lev <- nrow(real.simplex)
cong94 <- which(fulldb$congress==94)
real.simplex <- rbind(real.simplex,real.simplex[cong94,])
real.simplex[(lev+1):nrow(real.simplex),"group"] <- "94th Congress"

real.simplex <- real.simplex[complete.cases(real.simplex),]
half.c <- 1/sqrt(3)
real.simplex$Coord1 <- half.c*(real.simplex$Dissident-real.simplex$Partisan) 
pdf("../Graphs/ObsPorts.pdf",height=4.3,width=11.5)
xyplot(Loner~Coord1|group
       ,data=real.simplex
       ,layout=c(3,1)
       ,index.cond=list(c(3,1,2))
       ,xlim=c(-.58-0.25,.58+0.25)
       ,ylim=c(-0.15,1.15)
       ,xlab=NULL
       ,ylab=NULL
       ,half.c=1/sqrt(3)
       ,cex=c(1.5,0.7,0.7)
       ,col=c("black",rep(rgb(0,0,0,0.5),2))
       ,par.settings=standard.theme("pdf", color=FALSE)
       ,scales=list(x=list(at=NULL),y=list(at=NULL))
       ,panel=function(x,y,half.c,cex,col,...){
         if(panel.number()==1)
           panel.text(c(0),c(0.45),"Bill Young (R), 94th")
         panel.xyplot(x,y,pch=19,col=col[panel.number()],cex=cex[panel.number()])  
         panel.polygon(c(-half.c,0,half.c),c(0,1,0),lwd=1.5)
         panel.text(c(-half.c-.05,0,half.c+.05),c(-.05,1.05,-.05),c("Partisan","Loner","Dissident"))
       }
)
dev.off()



##Preprocess data
fulldb <- fulldb[order(fulldb$mc,fulldb$congress),]
fulldb$MeanDevFromMedian <- ave(fulldb$devFromMedian,fulldb$mc,FUN=function(x)mean(x,na.rm=TRUE))
fulldb$LagSimplexOpp <- ave(fulldb$AvgSimplexOpposition,fulldb$mc,FUN=function(x)c(NA,x[-length(x)]))
fulldb$LagSimplexCop  <-  ave(fulldb$AvgSimplexCopartisan,fulldb$mc,FUN=function(x)c(NA,x[-length(x)]))
fulldb$LagSimplexAlone  <-  ave(fulldb$AvgSimplexAlone,fulldb$mc,FUN=function(x)c(NA,x[-length(x)]))

fulldb.an <- subset(fulldb,select=c(TotCopartisans
                                    ,TotOpposition
                                    ,TotAlone
                                    ,nonAlone
                                    ,marvicprim
                                    ,marvicgen
                                    ,compStateSize
                                    ,MeanDevFromMedian
                                    ,devFromMedian
                                    ,mapercpres
                                    ,InMajority
                                    ,mc
                                    ,RedistYear
                                    ,leadership
                                    ,congress
                                    ,party
                                    ,marvicprimNext
                                    ,marvicgenNext
                                    ,seniority
                                    ,nSponsoredBills
                                    ,AvgSimplexCopartisan
                                    ,AvgSimplexOpposition
                                    ,AvgSimplexAlone
                                    ,Loner
                                    ,Partisan
                                    ,Dissident
                                    ,LagSimplexOpp
                                    ,LagSimplexCop
                                    ,LagSimplexAlone
))
fulldb.an <- subset(fulldb.an,complete.cases(fulldb.an[,-c(17,18,27:29)]))
fulldb.an <- fulldb.an[order(fulldb.an$congress),]
fulldb.an$post96 <- fulldb.an[,"congress"] > 96
fulldb.an$congress <- as.factor(fulldb.an$congress)
fulldb.an$party <- as.factor(fulldb.an$party)
fulldb.an$compStateSize  <- scale(fulldb.an$compStateSize)
fulldb.an$MeanDevFromMedian <- scale(fulldb.an$MeanDevFromMedian)
fulldb.an$devFromMedian <- scale(fulldb.an$devFromMedian)
fulldb.an$mapercpres <- scale(fulldb.an$mapercpres)
fulldb.an$LagSimplexOpp <- scale(fulldb.an$LagSimplexOpp)
fulldb.an$LagSimplexCop <- scale(fulldb.an$LagSimplexCop)
fulldb.an$LagSimplexAlone <- scale(fulldb.an$LagSimplexCop)
fulldb.an$seniority <- scale(fulldb.an$seniority)


#############################
## Restrospective Analysis
#############################

## Bayesian analysis of compositional data
## using mixed-effects KK multivariate t model.

#For debugging
#db_sample <- sample(1:nrow(cosp_simplex),100)
#fulldb.an <- fulldb.an[db_sample,]

des_mat <- model.matrix(~marvicprim*marvicgen 
                        +compStateSize
                        +devFromMedian
                        +mapercpres
                        +leadership
                        +RedistYear
                        +LagSimplexCop
                        +LagSimplexOpp
                        +seniority
                        +congress
                        ,data=fulldb.an)
fulldb.an_sub <- fulldb.an[rownames(fulldb.an)%in%rownames(des_mat),]
cosp_simplex <- with(fulldb.an_sub,cbind(Loner,Partisan,Dissident))+.001
cosp_simplex <- prop.table(cosp_simplex,1)

low_p_margin <- quantile(des_mat[,"marvicprim"],0)
low_g_margin <- quantile(des_mat[,"marvicgen"],0)
hi_p_margin <- quantile(des_mat[,"marvicprim"],1)
hi_g_margin <- quantile(des_mat[,"marvicgen"],1)

##Form prediction matrix
#For scenario predictions at mean values
#and for average predictive values
X_pred <- des_mat[1,]
for(i in 2:length(X_pred)) X_pred[i] <- 0 # Set predictors at mean
X_pred <- matrix(X_pred,nrow=4,ncol=length(X_pred),byrow = TRUE,dimnames=list(NULL,names(X_pred)))
cong <- c("congress106","congress102","congress107","congress95")
for(i in 1:4){
  mars <- switch(i
                 ,c(low_p_margin,low_g_margin)    
                 ,c(low_p_margin,hi_g_margin)
                 ,c(hi_p_margin,low_g_margin)
                 ,c(hi_p_margin,hi_g_margin)
  )
  X_pred[i,"marvicprim"] <- mars[1] 
  X_pred[i,"marvicgen"] <- mars[2]
  X_pred[i,"marvicprim:marvicgen"] <- mars[1]*mars[2]
  X_pred[i,cong[i]] <- 1
}
for(i in 1:4){
  mars <- switch(i
                 ,c(low_p_margin,low_g_margin)    
                 ,c(low_p_margin,hi_g_margin)
                 ,c(hi_p_margin,low_g_margin)
                 ,c(hi_p_margin,hi_g_margin)
  )
  temp_data <- transform(des_mat,"marvicprim"=mars[1],"marvicgen"=mars[2],"marvicprim.marvicgen"=mars[1]*mars[2])
  colnames(temp_data) <- colnames(X_pred)
  X_pred <- rbind(X_pred,temp_data)
}

stan_data <- list(N=nrow(des_mat)
                  ,N_pred = nrow(X_pred)
                  ,L=length(unique(fulldb.an_sub$mc))
                  ,P=ncol(des_mat)
                  ,T_full=ncol(cosp_simplex)
                  ,T=ncol(cosp_simplex)-1
                  ,Y=cosp_simplex
                  ,X=des_mat
                  ,leg=as.numeric(as.factor(fulldb.an_sub$mc))
                  ,X_pred = t(X_pred)
                  ,fresh=0
)


kk_model <- stan(file="Compositional.stan"
                 ,data=stan_data
                 ,chains=10
                 ,iter=1e3
                 ,cores=10)

save(kk_model,file="KKModel.RData")

load("KKModel.RData")

## Predictions of portfolios 

all_preds <- extract(kk_model,"Port_pred")
all_preds <- all_preds$Port_pred


## Under different scenarios
library(MASS)
scene=c("Both tough","Tough primary, easy general","Easy primary, tough general","Both easy")
scen_preds <- all_preds[,,1:4]
scen_pred_stack <- scen_preds[,,1]
for(i in 2:4) scen_pred_stack <- rbind(scen_pred_stack,scen_preds[,,i])
scen_pred_stack <- as.data.frame(scen_pred_stack)
scen_pred_stack$scene <- as.factor(rep(scene,each=nrow(scen_preds)))
scen_pred_stack$coord1 <- half.c*(scen_pred_stack[,3] - scen_pred_stack[,2]) 
scen_pred_stack$coord2 <- scen_pred_stack[,1] 
pdf("../Graphs/PredictedPorts.pdf",height=6,width=6.5)
xyplot(coord2~coord1|scene
       ,data=scen_pred_stack
       ,layout=c(2,2)
       ,index.cond=list(c(3,2,1,4))
       ,xlim=c(-.4-0.25,.3+0.25)
       ,ylim=c(-0.15,0.6062178+.15)
       ,xlab=NULL
       ,ylab=NULL
       ,par.settings=standard.theme("pdf", color=FALSE)
       ,scales=list(x=list(at=NULL),y=list(at=NULL))
       ,panel=function(x,y,...){
         z <-  kde2d(x,y)
         grid <- expand.grid(z$x,z$y)
         grid$z <- c(z$z)
         panel.polygon(c(-.4,-.05,.3),c(0,0.6062178,0),lwd=1.5)
         panel.contourplot(grid$Var1,grid$Var2,grid$z
                           ,subscripts=TRUE
                           ,contour=TRUE
                           ,region=FALSE
                           ,lwd=2.5
                           ,at=quantile(grid$z,probs=c(.95,.99))
                           ,...)
         panel.text(c(-.45,-0.05,.35),c(-.05,0.65,-.05),c("Partisan","Loner","Dissident"))
       }
)
dev.off()

## Distribution of predictive effects (as percent change)
first_end <- nrow(des_mat)+4
#From easy primary to tough primary, given easy general, on copartisan
pe_preds_cop <- 100*(all_preds[,,(first_end+1):(first_end+nrow(des_mat))] - all_preds[,,(first_end+nrow(des_mat)*2+1):(first_end+nrow(des_mat)*3)]
)/all_preds[,,(first_end+nrow(des_mat)*2+1):(first_end+nrow(des_mat)*3)]
pe_preds_cop <- colMeans(pe_preds_cop[,2,])

#From easy general to tough general, given easy primary, on opposition
pe_preds_opp <- 100*(all_preds[,,(first_end+nrow(des_mat)+1):(first_end+nrow(des_mat)*2)] - all_preds[,,(first_end+nrow(des_mat)*2+1):(first_end+nrow(des_mat)*3)]
)/all_preds[,,(first_end+nrow(des_mat)*2+1):(first_end+nrow(des_mat)*3)]
pe_preds_opp <- colMeans(pe_preds_opp[,3,])

#From two tough elections to two easy elections, on alone
pe_preds_al <- 10*(all_preds[,,5:first_end] - all_preds[,,(first_end+nrow(des_mat)*2+1):(first_end+nrow(des_mat)*3)]
)/all_preds[,,(first_end+nrow(des_mat)*2+1):(first_end+nrow(des_mat)*3)]
pe_preds_al <- colMeans(pe_preds_al[,1,])

all_fx <- data.frame(val=c(pe_preds_cop,pe_preds_opp,pe_preds_al)
                     ,grp=rep(c("Tougher primary,\ngiven easy general"
                                ,"Tougher general,\ngiven easy primary"
                                ,"Both get easier"),each=length(pe_preds_cop)))
levels(all_fx$grp) <- c("Both get easier","Tougher primary,\ngiven easy general","Tougher general,\ngiven easy primary")

pdf("../Graphs/FXDist.pdf",height=4.5,width=5)
boxplot(val ~ grp,data=all_fx,outline=FALSE
        ,main="Effects of Change in Electoral Challenge"
        ,ylab="Percent Change in Probability"
        ,cex.axis=0.7
        ,ylim=c(-10,60)
        ,frame=FALSE
        ,at=c(1,3,2)
) 
abline(h=0,lty=3)
text(c(1,3,2),c(52,59,37)
     ,c("Of working\nalone","Of working\nw. opposition","Of working\nw. copartisans")
     ,cex=0.7)
dev.off()
##Regression table for appendix
sd_alphas <- sd(extract(kk_model,"alpha_l")[[1]])
nu <- median(extract(kk_model,"nu")[[1]])
sigma_y <- extract(kk_model,"Sigma")[[1]]
all_coef <- round(summary(kk_model,c("beta"),c(0.025,0.5,0.975))$summary,3)[,c(4:7)]
reg_tab <- cbind(all_coef[1:stan_data$P,2]
                 ,aaply(all_coef[1:stan_data$P,c(1,3)],1,paste,collapse=",")
                 ,round(all_coef[1:stan_data$P,4])
                 ,all_coef[(stan_data$P+1):(nrow(all_coef)),2]
                 ,aaply(all_coef[(stan_data$P+1):(nrow(all_coef)),c(1,3)],1,paste,collapse=",")
                 ,round(all_coef[(stan_data$P+1):(nrow(all_coef)),4])
)
rownames(reg_tab) <- colnames(des_mat)
colnames(reg_tab) <- rep(c("Point Estimate","Credible Interval (95%)","Nr. effective samples"),2)
print(xtable(reg_tab),file="MainTable.tex")



########################
## Prospective analysis
########################

fulldb.prosp <- subset(fulldb.an
                       ,complete.cases(subset(fulldb.an
                                              ,select=c(marvicprim
                                                        ,marvicgen
                                                        ,Loner
                                                        ,Dissident
                                                        ,Partisan
                                                        ,mapercpres
                                                        ,compStateSize
                                                        ,MeanDevFromMedian
                                                        ,seniority
                                                        ,party
                                                        ,leadership
                                                        ,congress
                                                        ,marvicprimNext
                                                        ,marvicgenNext
                                              ))))
fulldb.prosp$party <- ifelse(fulldb.prosp$party==200,0,1)

eq1 <- marvicprimNext ~ Partisan*marvicprim*marvicgen+
  Dissident*marvicprim*marvicgen+
  compStateSize+
  mapercpres+
  MeanDevFromMedian+
  seniority+
  party+
  leadership+
  congress

eq2 <- update(eq1,marvicgenNext~.) 
post.model.prim <- systemfit(list(Primary=eq1,General=eq2)
                             ,method="SUR"
                             ,data=fulldb.prosp)
ci_pmp <- confint(post.model.prim)
cis <- paste("(",paste(round(ci_pmp[,1],3),round(ci_pmp[,2],3),sep=", "),")",sep="")
pes <- round(coef(post.model.prim),3)
print(xtable(cbind(pes[1:30],cis[1:30],pes[31:60],cis[31:60])) ,file="SUR.tex")

post_pred_x <- data.frame(marvicprim = rep(c(low_p_margin,hi_p_margin),each = 3)
                          ,marvicgen = rep(c(hi_g_margin,low_g_margin),each = 3)
                          ,Partisan =  c(.75, .27 ,0.45,.27,.75 ,.45)
                          ,Dissident = c(.25,0.46,0.45,.46,.25,.45)
                          ,compStateSize = mean(fulldb.prosp$compStateSize)
                          ,mapercpres = mean(fulldb.prosp$mapercpres)
                          ,MeanDevFromMedian = mean(fulldb.prosp$MeanDevFromMedian)
                          ,seniority = mean(fulldb.prosp$seniority)
                          ,party = 0
                          ,leadership = 0
                          ,congress = factor(101,levels=levels(fulldb.prosp$congress))
)
post_preds <- predict(post.model.prim,newdata=post_pred_x,interval="confidence",level=.9)
relevant_preds <- as.data.frame(rbind(as.matrix(post_preds[1:3,1:3]),as.matrix(post_preds[4:6,4:6])))
names(relevant_preds) <- c("Pred","LB","UB")
relevant_preds$scen <- rep(c("Tough primary, easy general (best: Partisan)","Easy primary, tough general (best: Dissident)"),each=3)
relevant_preds$port <- c("Dissident","Partisan","Loner","Dissident","Partisan","Loner")

pdf("../Graphs/ProspRes.pdf",height=3.5,width=8.5)
stripplot(port~Pred|scen
          ,data=relevant_preds
          ,UB=relevant_preds$UB
          ,LB=relevant_preds$LB
          ,index.cond=list(c(2,1))
          ,par.settings=standard.theme("pdf", color=FALSE)
          ,scales=list(x="free")
          ,xlab="Margin of victory in next election"
          ,xlim=list(c(0.01,.26),c(0.5,.95))
          ,panel=function(x,y,subscripts,LB,UB,...){
            panel.stripplot(x,y,pch=19,cex=1.2,...)
            panel.arrows(x0=LB[subscripts]
                         ,x1=UB[subscripts]
                         ,y0=y,y1=y
                         ,lwd=2,code=3,angle=90
                         ,length=0.05)
          })
dev.off()



#### Appendix analyses:

#Prospective Jacobson 2015 incumbency model
#augmented. 
fulldb <- fulldb[order(fulldb$state,fulldb$cd,fulldb$year),]
fulldb$DemShare <- with(fulldb,ifelse(party==100,marvicgen,-marvicgen))*100
#fulldb$DemShareLag <- with(fulldb,ave(DemShare,state,cd,FUN=function(x)c(NA,x[-length(x)])))
fulldb$partyGK <- ifelse(fulldb$party==100,1,-1)
fulldb$PresPercJ <- ifelse(fulldb$partyGK==1,fulldb$mapercpres,1-fulldb$mapercpres)*100
fulldb$PartisanLag <- with(fulldb,ave(Partisan,state,cd,FUN=function(x)c(NA,x[-length(x)])))
fulldb$DissidentLag <- with(fulldb,ave(Dissident,state,cd,FUN=function(x)c(NA,x[-length(x)])))
fulldb$incumbGK <- with(fulldb,ave(partyGK,state,cd,FUN=function(x)c(NA,x[-length(x)])))
fulldb$PartisanLag <- with(fulldb,ifelse(incumbGK==1,PartisanLag,-PartisanLag))
fulldb$DissidentLag <- with(fulldb,ifelse(incumbGK==1,DissidentLag,-DissidentLag))
j_model <- lm(DemShare~PresPercJ*(PartisanLag+DissidentLag)+partyGK+incumbGK+as.factor(congress),
              data=fulldb)
summary(j_model)
stargazer(j_model,out="jacobson.tex")

#Retrospective: models on freshmen only
fresh_mat <- model.matrix(~marvicprim*marvicgen 
                          +compStateSize
                          +mapercpres
                          +congress
                          ,data=subset(fulldb.an_sub,seniority==min(fulldb.an_sub$seniority)))
fresh_simplex <- with(subset(fulldb.an_sub,seniority==min(fulldb.an_sub$seniority)),cbind(Loner,Partisan,Dissident))+.001
fresh_simplex <- prop.table(fresh_simplex,1)
fresh_pred <- fresh_mat[1,]
for(i in 1:4){
  mars <- switch(i
                 ,c(low_p_margin,low_g_margin)    
                 ,c(low_p_margin,hi_g_margin)
                 ,c(hi_p_margin,low_g_margin)
                 ,c(hi_p_margin,hi_g_margin)
  )
  temp_data <- transform(fresh_mat,"marvicprim"=mars[1],"marvicgen"=mars[2],"marvicprim.marvicgen"=mars[1]*mars[2])
  colnames(temp_data) <- names(fresh_pred)
  fresh_pred <- rbind(fresh_pred,temp_data)
}
fresh_pred <- fresh_pred[-1,]

fresh_data <- list(N=nrow(fresh_mat)
                   ,N_pred = nrow(fresh_pred)
                   ,L=2
                   ,P=ncol(fresh_mat)
                   ,T_full=ncol(fresh_simplex)
                   ,T=ncol(fresh_simplex)-1
                   ,Y=fresh_simplex
                   ,X=fresh_mat
                   ,X_pred = t(fresh_pred)
                   ,leg=1:nrow(fresh_mat)
                   ,fresh=1
)

fresh_model <- stan(file="Compositional.stan"
                    ,data=fresh_data
                    ,chains=5
                    ,iter=7e2
                    ,cores=5)
save(fresh_model,file="FreshModel.RData")
load("FreshModel.RData")
all_preds_fresh <- extract(fresh_model,"Port_pred")
all_preds_fresh <- all_preds_fresh$Port_pred

## Distribution of predictive effects (as percent change)
first_end <- nrow(fresh_mat)
#From easy primary to tough primary, given easy general, on copartisan
pe_preds_cop_f <- 100*(all_preds_fresh[,,(first_end+1):(first_end+nrow(fresh_mat))] - all_preds_fresh[,,(first_end+nrow(fresh_mat)*2+1):(first_end+nrow(fresh_mat)*3)]
)/all_preds_fresh[,,(first_end+nrow(fresh_mat)*2+1):(first_end+nrow(fresh_mat)*3)]
pe_preds_cop_f <- colMeans(pe_preds_cop_f[,2,])

#From easy general to tough general, given easy primary, on opposition
pe_preds_opp_f <- 100*(all_preds_fresh[,,(first_end+nrow(fresh_mat)+1):(first_end+nrow(fresh_mat)*2)] - all_preds_fresh[,,(first_end+nrow(fresh_mat)*2+1):(first_end+nrow(fresh_mat)*3)]
)/all_preds_fresh[,,(first_end+nrow(fresh_mat)*2+1):(first_end+nrow(fresh_mat)*3)]
pe_preds_opp_f <- colMeans(pe_preds_opp_f[,3,])

#From two tough elections to two easy elections, on alone
pe_preds_al_f <- 10*(all_preds_fresh[,,1:first_end] - all_preds_fresh[,,(first_end+nrow(fresh_mat)*2+1):(first_end+nrow(fresh_mat)*3)]
)/all_preds_fresh[,,(first_end+nrow(fresh_mat)*2+1):(first_end+nrow(fresh_mat)*3)]
pe_preds_al_f <- colMeans(pe_preds_al_f[,1,])

all_fx_fresh <- data.frame(val=c(pe_preds_cop_f,pe_preds_opp_f,pe_preds_al_f)
                           ,grp=rep(c("Tougher primary,\ngiven easy general"
                                      ,"Tougher general,\ngiven easy primary"
                                      ,"Both get easier"),each=length(pe_preds_cop_f)))

pdf("../Graphs/FXDistFresh.pdf",height=4.5,width=5.5)
boxplot(val ~ grp,data=all_fx_fresh,outline=FALSE
        ,main="Effects of Change in Electoral Challenge (Freshmen)"
        ,ylab="Percent Change in Probability"
        ,cex.axis=0.7
        ,ylim=c(-10,60)
        ,frame=FALSE
        ,at=c(1,3,2)
) 
abline(h=0,lty=3)
text(c(1,3,2),c(17,52,33)
     ,c("Of working\nalone","Of working\nw. opposition","Of working\nw. copartisans")
     ,cex=0.7)
dev.off()

##Regression table for appendix
all_coef <- round(summary(kk_model,c("beta"),c(0.025,0.5,0.975))$summary,3)[,c(4:7)]
reg_tab <- cbind(all_coef[1:stan_data$P,2]
                 ,aaply(all_coef[1:stan_data$P,c(1,3)],1,paste,collapse=",")
                 ,round(all_coef[1:stan_data$P,4])
                 ,all_coef[(stan_data$P+1):(nrow(all_coef)),2]
                 ,aaply(all_coef[(stan_data$P+1):(nrow(all_coef)),c(1,3)],1,paste,collapse=",")
                 ,round(all_coef[(stan_data$P+1):(nrow(all_coef)),4])
)
rownames(reg_tab_fresh) <- colnames(des_mat)
colnames(reg_tab_fresh) <- rep(c("Point Estimate","Credible Interval (95%)","Nr. effective samples"),2)



#### PCA stuff

setwd("~/MyCloud/Cosponsorship/Rollcall/")

## Cosponsorship and Rollcalls using PCA
pca.ideal <-function(temp.data){
  Tt<-as.matrix(temp.data[,-c(1:7)])%*%t(as.matrix(temp.data[,-c(1:7)])) 
  Aff <- Tt
  Tt<-Tt/diag(Tt)
  Tt[Tt=="NaN"] <- 0
  diag(Tt) <- 1
  Tt<- sqrt(Tt)
  pMat <- Aff/sum(diag(Aff))
  PCTt <- prcomp(Tt,retx=TRUE)
  IdealTt <- PCTt$rotation[,1]
  IdealTt <- (2*((IdealTt-min(IdealTt))/(max(IdealTt)-min(IdealTt)))-1)
  color <- temp.data[,6]
  to.rotate<-lm(IdealTt~-1+factor(color))
  if (to.rotate$coefficients[2]<to.rotate$coefficients[1]) {IdealTt<-IdealTt*-1}
  party.medians <- by(IdealTt,color,median)
  dist.to.medians <- ifelse(color==10001|color==100,abs(IdealTt-party.medians[2]),abs(IdealTt-party.medians[1]))
  return(dist.to.medians)
}

col.vec <- c(3,5,2,2,8,5,11)
IdealPointDistances <- list()

for(i in 93:108){
  file <- paste("DIPUTADOS",i,".dta", sep="")
  data.cosp <- read.dta(file)
  data.cosp <- data.cosp[,-c(3,4,7:23)]
  data.cosp <- cbind(data.cosp[,c(4,2)],NA,NA,NA,data.cosp[,3],NA,data.cosp[,-c(1:4)])
  #Roll Call
  fileName <- paste("hou",i,"kh.ord",sep="")
  raw.data <- as.matrix(readLines(con=fileName))
  voteData <- do.call(rbind,lapply(
    apply(substring(raw.data, 37)
          ,1,strsplit,split=""),unlist)
  )
  voteData <- apply(voteData,2,recode,
                    recodes="1:3=1;else=0"
  )
  
  non.vote.u <- substr(raw.data, 0,36)
  non.vote <- do.call(rbind,apply(non.vote.u,
                                  1,
                                  function(x){
                                    el.data <- read.fwf(textConnection(x)
                                                        ,col.vec)
                                    closeAllConnections() 
                                    return(el.data)
                                  })) 
  temp.data.rc <- cbind(non.vote,voteData)
  temp.data.rc <- subset(temp.data.rc,temp.data.rc[,6]==10001|temp.data.rc[,6]==20001)
  ID.rc <- as.numeric(paste(temp.data.rc[,2],i,sep=""))
  ID.cosp <- as.numeric(paste(data.cosp[,2],i,sep=""))
  
  pca.ideal.rc <- pca.ideal(temp.data.rc)
  pca.ideal.rc <- cbind(ID.rc,pca.ideal.rc)
  colnames(pca.ideal.rc) <- c("MCCongressID","Distance.RC")
  pca.ideal.cosp <- pca.ideal(data.cosp)
  pca.ideal.cosp <- cbind(ID.cosp,pca.ideal.cosp)
  colnames(pca.ideal.cosp) <- c("MCCongressID","Distance.Cosp")
  pca.ideal.all <- merge(pca.ideal.rc,pca.ideal.cosp)
  
  
  IdealPointDistances[[i]] <- pca.ideal.all
}

IdealPointDistances <- do.call(rbind,IdealPointDistances)

fullrccosp <- merge(IdealPointDistances,fulldb)

cond.eff <- function(elModel){
  modVcov <- vcov(elModel)
  print(modVcov[13,13])
  varPrim <-modVcov[2,2]
  varGen <-modVcov[3,3]
  varInter <-modVcov[13,13]
  covInterPrim <- modVcov[13,2]
  covInterGen <- modVcov[13,3]
  maxGen <- quantile(elModel$model[,3],1,na.rm=TRUE)
  maxPrim <- quantile(elModel$model[,2],1,na.rm=TRUE)
  ef.maxgen.prim <- coef(elModel)[2] + coef(elModel)[13]*maxGen
  se.maxgen.prim <- sqrt(varPrim + ((maxGen)^2)*(varInter) + 2*maxGen*covInterPrim)
  ef.maxprim.gen <- coef(elModel)[3]+ coef(elModel)[13]*maxPrim
  se.maxprim.gen <- sqrt(varGen + ((maxPrim)^2)*(varInter) + 2*maxPrim*covInterGen)
  result <- cbind(rbind(ef.maxgen.prim,ef.maxprim.gen),rbind(se.maxgen.prim,se.maxprim.gen))
  rownames(result) <- c("PrimaryGetsTougher|EasyGeneral","GeneralGetsTougher|EasyPrimary")
  colnames(result) <- c("Effect","S.E.")
  return(result)
}

##Roll Call
rc.model <- lm(sqrt(Distance.RC)~I(1-marvicprim)*I(1-marvicgen) 
               +scale(compStateSize)
               +scale(mapercpres)
               +scale(leadership)
               +as.factor(congress)
               ,data=fullrccosp)
cond.eff(rc.model)

## Cosponsor Call
cp.model <- lm(sqrt(Distance.Cosp)~I(1-marvicprim)*I(1-marvicgen) 
               +scale(compStateSize)
               +scale(mapercpres)
               +scale(leadership)
               +as.factor(congress)
               ,data=fullrccosp)
cond.eff(cp.model)


