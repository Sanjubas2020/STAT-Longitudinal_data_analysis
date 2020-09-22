
/*qPCR data analysis in SAS using GLIMMIX*/
data crd;
input trt$  @@;
do obs = 1 to 4;
  input y @@;
  output;
end;
datalines;
Water     0.92 1.22 1.05 0.88
GFP       1.22 1.38 1.55 0.92
LacZ      1.18 1.19 1.43 0.58
COPE      1.06 1.7  1.08  1.36
;
proc print data=crd;
run;
proc glimmix data=crd;
class trt;
model y=trt / dist=gamma link=log;
lsmeans trt/diff lines cl plot=diffplot plot=meanplot(cl);lsmeans trt/ diff adjust=tukey;
*lsmeans trt/diff=control('gfp') adjust=dunnett;
run;
/*proportion data best fit the beta distribution*/
/*log provides positive value*/

/*SAS codes for vATPase survival data*/
/* Data is presented in Multivariate form                 */
data multivar;
 input trt subj t1-t6; 
datalines;
1 1 100  95.8	75.0	62.5	58.3	54.2
1 2 100  100    91.7	83.3	79.2	66.7
1 3 100  83.3	79.2	62.5	50.0	41.7
2 1 100  83.3   75.0    62.5	54.2	41.7 
2 2 100  95.8	87.5	87.5	83.3	79.2
2 3 100  91.7	79.2	75.0	62.5	50.0
3 1 100  91.7	79.2	62.5	54.2	41.7
3 2 100  91.7	75.0	66.7	54.2	45.8
3 3 100  83.3	83.3	75.0	66.7	66.7

/* convert multivariate to univariate form, required for Proc MIXED & GLIMMIX        */
data univariate;
 set multivar;
 time=1; y=t1; output;
 time=2; y=t2; output;
 time=3; y=t3; output;
 time=4; y=t4; output;
 time=5; y=t5; output;
 time=6; y=t6; output;

proc print data=univariate; run;
proc print data=multivar; run;


/Determine what correlation structure do you have?? */
/* AR(1) ? */
title AR1;
proc glimmix  data=univariate;
  class subj trt time;
  model y=trt|time; 
  random intercept /subject=subj(trt);
  random time / subject=subj(trt) type=ar(1) residual;
 run;
 title;
/* ANTE(1) ? */
 title Ante;
proc glimmix  data=univariate;
  class subj trt time;
  model y=trt|time; 
  random time / subject=subj(trt) type=ante(1) residual;
 run;
 title; 
/* UN ? */
  title un;
proc glimmix  data=univariate;
  class subj trt time;
  model y=trt|time; 
  random time / subject=subj(trt) type=un residual;
 run;
 title;

/* CS ? */
 title CS;
proc glimmix  data=univariate;
  class subj trt time;
  model y=trt|time; 
  random time / subject=subj(trt) type=cs residual;
 run;
 title;
/* AR(1) wins */
ods graphics on;
 proc glimmix  data=univariate;
  class subj trt time;
  model y= trt trt|time/ddfm=kr; 
  random intercept/subject=subj(trt);
  random time / subject=subj(trt) type=ar(1) residual;
  lsmeans trt*time/plot=meanplot(sliceby=trt join);
/*  lsmeans program/diff=control('CONT') adjust=dunnett;*/
  lsmeans time/diff adjust=tukey;
 run;
