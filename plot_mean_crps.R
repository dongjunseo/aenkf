########################################################################
# plot_mean_crps.R - version Aug 14, 2021, by D.-J. Seo at UTA/HWRL
########################################################################
library(hydroGOF)

X11(width=10,height=7,pointsize=14)

x<-0;crps_e<-0;crps_a<-0;crps_b<-0;crps_aenkf<-0;crps_s<-0
rmse_sim<-0;rmse_DA<-0;rmse_ens<-0

for(kk in  1:4) {
if(kk == 1) {cid<-'gtbm3';cjd<-'GTBM3';rfc='NE'}
if(kk == 2) {cid<-'coli2';cjd<-'COLI2';rfc='NC'}
if(kk == 3) {cid<-'dltc1';cjd<-'DLTC1';rfc='CN'}
if(kk == 4) {cid<-'monn7';cjd<-'MONN7';rfc='SE'}

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

outfile1<-paste('mean_crps_',cjd,'.jpg',sep='')
infile2<-paste('./output',id,'/',cjd,'/run_',cjd,'_aenkf_',ci,cj,'/',cjd,'.fcst_suppl',sep='')
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
ctitle<-paste(rfc,' - ',cjd,sep='')

infile1<-paste('./output',id,'/',cjd,'/run_',cjd,'_enkf_',ci,cj,'/',cjd,'.fcst_obs',sep='')
infile2<-paste('./output',id,'/',cjd,'/run_',cjd,'_enkf_',ci,cj,'/',cjd,'.fcst_base',sep='')
infile3<-paste('./output',id,'/',cjd,'/run_',cjd,'_enkf_',ci,cj,'/',cjd,'.fcst_enkf',sep='')
infile4<-paste('./output',id,'/',cjd,'/run_',cjd,'_aenkf_',ci,cj,'/',cjd,'.fcst_aenkf',sep='')
infile5<-paste('./output',id,'/',cjd,'/run_',cjd,'_enkf_',ci,cj,'/',cjd,'.crps_base',sep='')
infile6<-paste('./output',id,'/',cjd,'/run_',cjd,'_enkf_',ci,cj,'/',cjd,'.crps_enkf',sep='')
infile7<-paste('./output',id,'/',cjd,'/run_',cjd,'_aenkf_',ci,cj,'/',cjd,'.crps_aenkf',sep='')
infile8<-paste('./output',id,'/',cjd,'/run_',cjd,'_enkf_',ci,cj,'/',cjd,'.fcst_sngl',sep='')

qobs<-read.table(infile1)
qsim<-read.table(infile2)
qenkf<-read.table(infile3)
qens<-qenkf
qcbpk<-read.table(infile4)
qDA<-qcbpk
crps_base<-read.table(infile5)
crps_enkf<-read.table(infile6)
crps_aenkf<-read.table(infile7)
crps_sv<-read.table(infile8)

size<-min(nrow(qobs),nrow(qsim),nrow(qenkf),nrow(qcbpk),nrow(qDA),nrow(qens),nrow(crps_base),nrow(crps_enkf),nrow(crps_aenkf))
 num_data<-dim(qDA)
#
# condition on rising vs. falling limb
#
 icut<-(-2)  #include all
#icut<-( 1)  #include rising limb only
ii<-icut
for (i in 1:num_data[2]){
 x[i]<-(i-1)*0.25
if(icut == -2) {
  qs1<-qsim[,i][irise>ii]
  qo1<-qobs[,i][irise>ii]
  qe1<-qens[,i][irise>ii]
  qd1<-qDA[,i][irise>ii]
  aopt1<-aopt[irise>ii]
  dfs1<-dfs[irise>ii]
  cb1<-crps_base[,i][irise>ii]
  ce1<-crps_enkf[,i][irise>ii]
  ca1<-crps_aenkf[,i][irise>ii]
  cs1<-crps_sv[,i][irise>ii]
}
if(icut != -2) {
  qs1<-qsim[,i][irise==ii]
  qo1<-qobs[,i][irise==ii]
  qe1<-qens[,i][irise==ii]
  qd1<-qDA[,i][irise==ii]
  aopt1<-aopt[irise==ii]
  dfs1<-dfs[irise==ii]
  cb1<-crps_base[,i][irise=ii]
  ce1<-crps_enkf[,i][irise=ii]
  ca1<-crps_aenkf[,i][irise=ii]
  cs1<-crps_sv[,i][irise=ii]
}
lqs1<-length(qs1[!is.na(qs1)])
lqo1<-length(qo1[!is.na(qo1)])
lqe1<-length(qe1[!is.na(qe1)])
lqd1<-length(qd1[!is.na(qd1)])
laopt1<-length(aopt1[!is.na(aopt1)])
ldfs1<-length(dfs1[!is.na(dfs1)])
lcb1<-length(cb1[!is.na(cb1)])
lce1<-length(ce1[!is.na(ce1)])
lca1<-length(ca1[!is.na(ca1)])
lcs1<-length(cs1[!is.na(cs1)])

size<-len<-min(lqs1,lqo1,lqe1,lqd1,lcb1,lce1,lca1,lcs1)
qs1<-qs1[1:size]
qo1<-qo1[1:size]
qe1<-qe1[1:size]
qd1<-qd1[1:size]
cb1<-cb1[1:size]
ce1<-ce1[1:size]
ca1<-ca1[1:size]
cs1<-cs1[1:size]
hsize<-1
ao1<-aopt1[hsize:size]
df1<-dfs1[hsize:size]

 ibeg<-1
 iend<-size
 qs1<-qs1[ibeg:iend]
 qo1<-qo1[ibeg:iend]
 qe1<-qe1[ibeg:iend]
 qd1<-qd1[ibeg:iend]
 ao1<-ao1[ibeg:iend]
 df1<-df1[ibeg:iend]
 cb1<-cb1[ibeg:iend]
 ce1<-ce1[ibeg:iend]
 ca1<-ca1[ibeg:iend]
 cs1<-cs1[ibeg:iend]
#
# condition on magnitude of alpha
#
 acut<-(-0.5   ) #include all
#acut<-( 0.5   ) #include only if alpha > 0.5
qs1<-qs1[ao1>acut];qe1<-qe1[ao1>acut];qd1<-qd1[ao1>acut];qo1<-qo1[ao1>acut];df1<-df1[ao1>acut];cb1<-cb1[ao1>acut];ce1<-ce1[ao1>acut];ca1<-ca1[ao1>acut];cs1<-cs1[ao1>acut]
##################################
# specify the cutoff
##################################
#qcut<-200
#if(kk == 6) qcut<-30
#if(kk == 9) qcut<-30
 qcut<-0 #include all flows
 qs1<-qs1[qo1>qcut]
 qe1<-qe1[qo1>qcut]
 qd1<-qd1[qo1>qcut]
 qo1<-qo1[qo1>qcut]
 cb1<-cb1[qo1>qcut]
 ce1<-ce1[qo1>qcut]
 ca1<-ca1[qo1>qcut]
 cs1<-cs1[qo1>qcut]

 dcut<-(1000000.)
 qs1<-qs1[df1<dcut]
 qe1<-qe1[df1<dcut]
 qd1<-qd1[df1<dcut]
 qo1<-qo1[df1<dcut]
 cb1<-cb1[df1<dcut]
 ce1<-ce1[df1<dcut]
 ca1<-ca1[df1<dcut]
 cs1<-cs1[df1<dcut]

 rmse_sim[i]<-rmse(qs1,qo1)
 rmse_DA[i]<-rmse(qd1,qo1)
 rmse_ens[i]<-rmse(qe1,qo1)

 crps_s[i]<-mean(abs(cs1-qo1))
 crps_b[i]<-mean(abs(cb1))
 crps_e[i]<-mean(ce1)
 crps_a[i]<-mean(ca1)
}

if(jj==0) {
par(mfrow=c(1,1))
par(pin=c(4.,4.))
par(lwd=1)
par(cex=1.)
xrange<-range(0,ix*0.25)
  ymin<-0
  ymax<-max(crps_b,crps_e,crps_a,crps_s)
  ymax<-ymax*1.0
  yrange<-range(ymin,ymax)
}
  if(jj==0) {
  plot(x[1:21],crps_s,xlim=xrange,ylim=yrange,xlab='',ylab='',type='l',lty=1,col=5,log='')
  par(new=TRUE)
  plot(x[1:21],crps_b,xlim=xrange,ylim=yrange,xlab='',ylab='',type='l',lty=1,col=3,log='')
  par(new=TRUE)
  plot(x[1:21],crps_a,xlim=xrange,ylim=yrange,xlab='',ylab='',type='l',lty=1,col=2,log='')
  par(new=TRUE)
  plot(x[1:21],crps_e,xlim=xrange,ylim=yrange,xlab='LEAD TIME (DAYS)',ylab='MEAN CRPS (CMS)',type='l',lty=1,col=4,log='')
redu_enkf0<-100.*(crps_b[1]-crps_e[1])/crps_b[1]
redu_DA0<-100.*(crps_b[1]-crps_a[1])/crps_b[1]
redu_enkf1<-100.*(crps_b[5]-crps_e[5])/crps_b[5]
redu_DA1<-100.*(crps_b[5]-crps_a[5])/crps_b[5]
print(c(redu_enkf0,redu_DA0,redu_enkf1,redu_DA1))
}
}
par(cex=0.90)
legend(0.5*(xmax-xmin)+xmin,0.4*(ymax-ymin)+ymin,legend=c("SV W/O DA","ENS W/O DA"," ","WC EnKF","WC AEnKF"),pch=c(" "," "," "," "," "," "," "," "),col=c(5,3,1,4,2,1,4,2,1),lty=c(1,1,0,1,1,0,2,2,0),bty="n",lwd=1.5)
title(ctitle)
 dev2bitmap(file=outfile1,type="jpeg",height=6,width=9,res=300,pointsize=12)
}

q()
