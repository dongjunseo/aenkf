#######################################################################
# plot_rmse.R - Version Aug 12, 2021, by D.-J. Seo at UTA/HWRL 
#######################################################################
X11(width=10,height=7,pointsize=14)

library(hydroGOF)

x<-0
for(kk in  1:4) {

if(kk == 1) {bnamel<-'gtbm3';bnameu<-'GTBM3';rfc='NE'}
if(kk == 2) {bnamel<-'coli2';bnameu<-'COLI2';rfc='NC'}
if(kk == 3) {bnamel<-'dltc1';bnameu<-'DLTC1';rfc='CN'}
if(kk == 4) {bnamel<-'monn7';bnameu<-'MONN7';rfc='SE'}

istep<-6
if(istep == 1) ix<-120
if(istep == 6) ix<-20

xmin<-0
xmax<-ix*0.25

ci<-'1'

for(jj in 0:0) {

if(jj == 0) {
cj<-'0'
id<-'0'
}

outfile1<-paste('rmse_',bnameu,'.jpg',sep='')

infile2<-paste('./output',id,'/',bnameu,'/run_',bnameu,'_aenkf_',ci,cj,'/',bnameu,'.fcst_suppl',sep='')
big<-scan(infile2,list(it=0,irise=0,ia=0,aopt=0,crps=0,qo=0,dfs=0,qda=0,q1=0,q10=0,q90=0,q99=0,qbe=0,nneg1=0,nneg2=0,drise=0,b1=0,b10=0,b90=0,b99=0,sumu_t=0,sumu1t=0,sumu_s=0))
it<-big$it
irise<-big$irise
ia<-big$ia
aopt<-big$aopt
crps<-big$crps
qo<-big$qo
dfs<-big$dfs
qda<-big$qda
qbe<-big$qbe
drise<-big$drise
ctitle<-paste(rfc,' - ',bnameu,sep='')
infile1<-paste('./output',id,'/',bnameu,'/run_',bnameu,'_enkf_',ci,cj,'/',bnameu,'.fcst_obs',sep='')
infile2<-paste('./output',id,'/',bnameu,'/run_',bnameu,'_enkf_',ci,cj,'/',bnameu,'.fcst_sngl',sep='')
infile3<-paste('./output',id,'/',bnameu,'/run_',bnameu,'_enkf_',ci,cj,'/',bnameu,'.fcst_enkf',sep='')
infile4<-paste('./output',id,'/',bnameu,'/run_',bnameu,'_aenkf_',ci,cj,'/',bnameu,'.fcst_aenkf',sep='')
infile5<-paste('./output',id,'/',bnameu,'/run_',bnameu,'_enkf_',ci,cj,'/',bnameu,'.fcst_base',sep='')
qobs<-read.table(infile1)
qsim<-read.table(infile2)
qenkf<-read.table(infile3)
qens<-qenkf
qcbpk<-read.table(infile4)
qDA<-qcbpk
qbase<-read.table(infile5)

size<-min(nrow(qobs),nrow(qsim),nrow(qenkf),nrow(qcbpk),nrow(qDA),nrow(qens),nrow(qbase))
#size<-(size-1)
 num_data<-dim(qDA)
 rmse_sim<-0;rmse_DA<-0;rmse_ens<-0;rmse_base<-0
#
# condition on rising vs. falling limb
#
 icut<-(-2) #include all
#icut<-( 1) #include rising limb only
ii<-icut
for (i in 1:num_data[2]){
x[i]<-(i-1)*0.25
if(icut == -2) {
 qs1<-qsim[,i][irise>ii]
 qo1<-qobs[,i][irise>ii]
 qe1<-qens[,i][irise>ii]
 qd1<-qDA[,i][irise>ii]
 qb1<-qbase[,i][irise>ii]
 aopt1<-aopt[irise>ii]
 dfs1<-dfs[irise>ii]
}
if(icut != -2) {
 qs1<-qsim[,i][irise==ii]
 qo1<-qobs[,i][irise==ii]
 qe1<-qens[,i][irise==ii]
 qd1<-qDA[,i][irise==ii]
 qb1<-qbase[,i][irise==ii]
 aopt1<-aopt[irise==ii]
 dfs1<-dfs[irise==ii]
}
lqs1<-length(qs1[!is.na(qs1)])
lqo1<-length(qo1[!is.na(qo1)])
lqe1<-length(qe1[!is.na(qe1)])
lqd1<-length(qd1[!is.na(qd1)])
lqb1<-length(qb1[!is.na(qb1)])
laopt1<-length(aopt1[!is.na(aopt1)])
ldfs1<-length(dfs1[!is.na(dfs1)])
size<-len<-min(lqs1,lqo1,lqe1,lqd1,lqb1)
qs1<-qs1[1:size]
qo1<-qo1[1:size]
qe1<-qe1[1:size]
qd1<-qd1[1:size]
qb1<-qb1[1:size]
hsize<-1
ao1<-aopt1[hsize:size]
df1<-dfs1[hsize:size]

 ibeg<-1
#ibeg<-2750

 iend<-size
 qs1<-qs1[ibeg:iend]
 qo1<-qo1[ibeg:iend]
 qe1<-qe1[ibeg:iend]
 qd1<-qd1[ibeg:iend]
 qb1<-qb1[ibeg:iend]
 ao1<-ao1[ibeg:iend]
 df1<-df1[ibeg:iend]
#
# condition on alpha
#
 acut<-(-0.5   ) #include all 
#acut<-( 0.5   ) #include only if alpha > 0.5
 qs1<-qs1[ao1>acut];qe1<-qe1[ao1>acut];qd1<-qd1[ao1>acut];qb1<-qb1[ao1>acut];qo1<-qo1[ao1>acut];df1<-df1[ao1>acut]
#qs1<-qs1[ao1==acut];qe1<-qe1[ao1==acut];qd1<-qd1[ao1==acut];qo1<-qo1[ao1==acut];df1<-df1[ao1==acut]

#qcut<-200
#if(kk == 6) qcut<-30
#if(kk == 9) qcut<-30
 qcut<-0
 qs1<-qs1[qo1>qcut]
 qe1<-qe1[qo1>qcut]
 qd1<-qd1[qo1>qcut]
 qb1<-qb1[qo1>qcut]
 qo1<-qo1[qo1>qcut]

 dcut<-(1000000.)
 qs1<-qs1[df1<dcut]
 qe1<-qe1[df1<dcut]
 qd1<-qd1[df1<dcut]
 qb1<-qb1[df1<dcut]
 qo1<-qo1[df1<dcut]

 rmse_sim[i]<-rmse(qs1,qo1)
 rmse_DA[i]<-rmse(qd1,qo1)
 rmse_base[i]<-rmse(qb1,qo1)
 rmse_ens[i]<-rmse(qe1,qo1)
}

if(jj==0) {
par(mfrow=c(1,1))
par(pin=c(4.00,4.00))
par(lwd=1)
par(cex=1.)
xrange<-range(0,ix*0.25)
ymin<-0
ymax<-max(rmse_sim,rmse_DA,rmse_ens)
ymax<-ymax*1.0
yrange<-range(ymin,ymax)
}
if(jj==0) {
plot(x[1:21],rmse_sim,xlim=xrange,ylim=yrange,xlab='',ylab='',type='l',lty=1,col=3,log='')
par(new=TRUE)
plot(x[1:21],rmse_base,xlim=xrange,ylim=yrange,xlab='',ylab='',type='l',lty=1,col=5,log='')
par(new=TRUE)
plot(x[1:21],rmse_DA,xlim=xrange,ylim=yrange,xlab='',ylab='',type='l',lty=1,col=2,log='')
par(new=TRUE)
plot(x[1:21],rmse_ens,xlim=xrange,ylim=yrange,xlab='LEAD TIME (DAYS)',ylab='RMSE (CMS)',type='l',lty=1,col=4,log='')
redu_enkf0<-100.*(rmse_base[1]-rmse_ens[1])/rmse_base[1]
redu_DA0<-100.*(rmse_base[1]-rmse_DA[1])/rmse_base[1]
redu_enkf1<-100.*(rmse_base[5]-rmse_ens[5])/rmse_base[5]
redu_DA1<-100.*(rmse_base[5]-rmse_DA[5])/rmse_base[5]
print(c(redu_enkf0,redu_DA0,redu_enkf1,redu_DA1))
}
}
par(cex=0.90)
legend(0.5*(xmax-xmin)+xmin,0.4*(ymax-ymin)+ymin,legend=c("SV W/O DA","Ens W/O DA"," ","WC EnKF","WC AEnKF"),pch=c(" "," "," "," "," "," "," "," "),col=c(5,3,1,4,2,1,4,2,1),lty=c(1,1,0,1,1,0,2,2,0),bty="n",lwd=1.5)
title(ctitle)

dev2bitmap(file=outfile1,type="jpeg",height=6,width=9,res=300,pointsize=12)

}

q()
