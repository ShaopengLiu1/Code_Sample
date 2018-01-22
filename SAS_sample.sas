# Shaopeng's SAS sample
# part1 (line 7-140) : SAS class final project: data cleaning, and modeling for SPB and DBP 
# part2 (line 150-244): SAS macro for bootstrap resampling
# part3 (line 250-319): SAS macro for recursive selection for optimization
# part4 (line 320-600): frequently used SAS procedures and command
###########################################
# Part1, SAS class final project:  
# 1.1, import file
data DEMOG;
infile "/home/gems/m21-503/final/demog.csv" dlm=',' firstobs=2 missover pad dsd;
length ID $12;
length race $32 studyno $16;			
informat DOB mmddyy10. visitdat mmddyy10.;
input gender $ dob studyno$ race $ visitdat age id$;
format DOB mmddyy10. visitdat mmddyy10.;
label race="race recorded by country, balck/white and hispanic" dob="date of birth";
run;

proc print data=demog(obs=20) label;
proc contents data=demog;
run;

# 1.2, merge correlated data
libname final "/home/gems/m21-503/final";

proc sort data=demog out=demog2;
by ID;
proc sort data=final.bp out=bp;
by id;
proc sort data=final.studydata out=studydata;
by id;
run;

data mydata2;
merge demog2(IN=in1) bp(IN=in2) studydata(IN=in3);
by id;
if in1 & in2 & in3;
run;

proc print data=mydata2(obs=20);
proc contents data=mydata2;
run;

# 1.3, data step edit
data mydata3;
set mydata2;
format _bmi 5.1 _map 6.1;
if age >= 50 then old=1;
else if 0 <= age <50 then old=2;
else old=.;
label old="50 or older as 1, less than 50 as 2";
whratio = waist / hip;
label whratio = "ratio of waist over hip";
drop waist hip;
run;

proc print data=mydata3(obs=20);
proc contents data=mydata3;
run;

# 1.4, add format
proc format;
value myfmt
1 = "Older Subjects"
2 = "Younger Subjects";
run;

data mydata5;
set mydata3;
format old myfmt.;
run;

# 1.5, data tarnsformation
data mydata6;
set mydata5;
if race="Black: Non Hispanic" or race="White: Non Hispanic";
rename _sbp=sbp _dbp=dbp;
new_race = substr(race, 1,5);
sum_alcohol = sum(beer, wine, liqr, cooler, sake);
run;

# 1.6, basic stats and plots with ODS
ODS rtf file="/home/shaopeng.liu/final9.rtf";

proc sort data=mydata5 out=mydata7;
by race;
run;
proc means data=mydata7;
var _sbp;
by race;
run;
/* The black subjects have higher average systolic blood pressure. */

proc sort data=mydata5 out=mydata8;
by studyno;
proc freq data=mydata8;
tables smokenow /missing;
by studyno;
run;
/* So the GenNet has the highest rate of smoking. */

ODS rtf close;

proc sgplot data=mydata6;
histogram whratio;
run;

* This result is for Black/White only;
ODS graphics on / imagename="Alcohol" imagefmt=jpg;
proc sgplot data=mydata6;
hbar sum_alcohol /missing;
run;
ODS graphics on / reset=ALL;

# 1.7, use macro to edit data
%macro calmean(data= , var=, out= );	*just store the result for convenience;
proc means data=&data N NMISS mean;
var &var;
output out=&out;
run;
%mend calmean;

* run by self;
%calmean(data=lsp.final, var=sbp, out=sbpout);
%calmean(data=lsp.final, var=dbp, out=dbpout);

proc print data=sbpout;
proc print data=dbpout;
run;

# 1.8, find correlation and model the variable
proc corr data=lsp.final;
var sbp dbp;
run;

proc reg data=lsp.final;
model sbp=dbp;
run;

###########################################







# Part2, bootstrap resampling
# Project: Bootstrap Resampling
%macro resample(k0,total=100,itera=2000,boots=100);

%include "/home/shaopengliu/Outlier_Detecting/MLE_Single.sas";

data GaSam(drop=N r m j k);
  call streaminit(123);
  N=&total;
  r=5;
  m=18;
  array iteration[1:%eval(&itera)];
  do id=1 to N;
    if id<=%eval(&total-&k0) then do j=1 to %eval(&itera);
      iteration[j]=2*rand('gamma',r);
    end;
    else do j=1 to %eval(&itera);
      iteration[j]=2*rand('gamma',2);
      do k=1 to (m-1);
        iteration[j]=iteration[j]+2*rand('gamma',2);
      end;
    end;
    output;
  end;
run;

data output; label counts1="Number of outliers 1" counts2="Number of outliers 2"
                   scale1="Estimated scale1" scale2="Estimated scale2";
run;

%do i=1 %to %eval(&itera);
  data iteration; set gasam; keep iteration%eval(&i);
  run;
  %let r1=0;
  
  proc sort data=gasam(keep=iteration%eval(&i)) out=iteration2; by iteration%eval(&i);
  data iteration2; set iteration2(obs=%eval(&total-10));run;
  %let r2=0;
  
  %do j=1 %to %eval(&boots);
    proc surveyselect data=iteration out=temp1 NOPRINT
      method=urs sampsize=%eval(&total-10) outhits;
    run;
    
    proc means noprint data=temp1;
      var iteration%eval(&i);
      output out=temp n=__n var=__var mean=__mean;
    run;
    
    data _null_;set temp; call symput('scale1',__n*__mean*__mean/__var/(__n-1));
    run;
    
    %let r1=%sysevalf(&r1+&scale1);
    
    proc surveyselect data=iteration2 out=temp2 NOPRINT
      method=urs sampsize=%eval(&total-10) outhits;
    run;
    
    proc means noprint data=temp2;
      var iteration%eval(&i);
      output out=temp n=__n var=__var mean=__mean;
    run;
    
    data _null_;set temp; call symput('scale2',__n*__mean*__mean/__var/(__n-1));
    run;
    
    %let r2=%sysevalf(&r2+&scale2);
  %end;
  
  %let r1=%sysevalf(&r1/&boots);
  %let r2=%sysevalf(&r2/&boots);
  
  data temp1; label counts="Proportion";
  data temp2; label counts="Proportion";
  run;
  
  %mle(iteration%eval(&i),k=10,alpha=0.05,r=&r1,data=gasam,data2=temp1,id=id);
  %mle(iteration%eval(&i),k=10,alpha=0.05,r=&r2,data=gasam,data2=temp2,id=id);
  
  data _null_; set temp1; call symput('k1',counts);
  data _null_; set temp2; call symput('k2',counts);
  run;

  data temp; counts1=&k1; counts2=&k2; scale1=&r1; scale2=&r2;
  data output; set output temp;
  run;
%end;

title "#outliers=&k0";
proc freq data=output; tables counts1 counts2;
proc means data=output; var scale1 scale2;
run;

%mend;
####################################################



# Part3, recursive selection for optimization
# Project: Recursive optimization
%macro recursive(k0,total=100,itera=2000);

%include "/home/shaopengliu/Outlier_Detecting/MLE_Single.sas";

data GaSam(drop=N r m j k);
  call streaminit(123);
  N=&total;
  r=5;
  m=18;
  array iteration[1:%eval(&itera)];
  do id=1 to N;
    if id<=%eval(&total-&k0) then do j=1 to %eval(&itera);
      iteration[j]=2*rand('gamma',r);
    end;
    else do j=1 to %eval(&itera);
      iteration[j]=2*rand('gamma',2);
      do k=1 to (m-1);
        iteration[j]=iteration[j]+2*rand('gamma',2);
      end;
    end;
    output;
  end;
run;

data output; label counts="Number of outliers" scale="Estimated scale";
run;

%do i=1 %to %eval(&itera);
  %let k=-1;
  %let scale=0;
  
  %do %until(&k>&ka);
    %let k=%eval(&k+1);
    %let scale0=&scale;
    proc sort data=gasam(keep=iteration%eval(&i)) out=removed ; by iteration%eval(&i);
    data removed; set removed(obs=%eval(&total-&k));run; 
      
    proc means noprint data=removed;
      var iteration%eval(&i);
      output out=temp n=__n var=__var mean=__mean;
    run;
  
    data _null_;set temp; call symput('scale',__n*__mean*__mean/__var/(__n-1));
    run;
    
    data temp1; label counts="Proportion";
    run;
    
    %mle(iteration%eval(&i),k=10,alpha=0.05,r=&scale,data=gasam,data2=temp1,id=id);
  
    data _null_; set temp1; call symput('ka',counts);
    run;
  %end;
  
  data temp2; counts=&ka; scale=&scale0;
  data output; set output temp2;
  run;
%end;

title "#outliers=&k0";
proc freq data=output; tables counts;
proc means data=output; var scale;
run;

%mend;


###################################################


# 1, overview:
# 1.1, Sample mean comparison: t-test and Wilcoxon's test
data one;
set two;
sbpdiff=sbp2-sbp1;
label sbpdiff='change in sbp from pre to post';

proc ttest; 
classes group;
var platelet;
title "Unpaired t-test on platelet data"

proc means n mean std stderr t prt;
var sbpdiff;
title "Paired t-test on sbp data";

proc npar1way wilcoxon; 
classes group;
var isi;
title "Wilcoxon test on infarct size data";
rn;

# 2, categorical analysis:
# 2.1, Binomial distribution and proportion estimate
proc freq data=<my_data>;
table var / binomial (level=I);  # use binomialc option for continuity correction
table var / binomial (p=0.5 level=I);  # compare proportion
table var / binomial (equiv p=0.5 level=I margin=0.1); #equivelent test
table var / binomial (sup p=0.5 level=I margin=0.1); #superiority test
table var / binomial (noninf p=0.5 level=I margin=0.1); #noninferiority example
title I "Proportion with complete response"
run;

# 2.2, Chi-square:
proc freq data=<my_data>;
table prefer / nocum chisq testp=(20 80); title "Chi-square goodness-of-fit test";
table prefer*flavor / nocum chisq; title "Chi-square test of Independence";
run;

# 2.3, Fisher Exact test:
proc freq data=<my_data>;
table flavor*prefer / nocum fisher;
title "Fisher's Exact test";
run;

# 2.4, ordinal category: trend
proc freq data=<my_data>;
table flavor * prefer / nocum trend jt;
title "CA and JT tests";
run;

# 2.5, Gamma and Kendall's tau
proc freq data=one;
table bp*weight /chisq measures cl;
title "gamma and kendalls tau for short ordinal scales";
run;

# 2.6, CMH: common risk and odds ratio
proc freq data=one;
tables gender * treatment * response / cmh;
title "Effect of Treatment on XXX controlling for gender";
run;

# 3, Linear model: proc reg and proc glm
# 3.1, ANOVA with 
proc format;
value marrige
1 = 'married'
2 = 'widowed'
3 = 'separated'
4 = 'divorced'
5 = 'never married';
run;

data one;
set oasis.baseline;
IF IDPROJI < 9000; * delete duplicate IDs;
dep_sx = (DLEAT=5) + (DLLW=5) + (DLGW=5) + (DLSLP=5) + (DLTMS=5) + (DLTIRE=5) + (DLSLOW=5) + (DLHYP=5) 
+ (DLSEX=5) + (DLGUILT=5) + (DLCONCEN=5) + (DLMIXED=5) + (DLDEAD=5) + (DLDIE=5) + (DLSUIC=5) + (DLATTP =5);
label dep_sx='# depression symptoms';
format marital marrige.;
keep dep_sx marital haspets age;
run;

proc glm data=one order=formatted; 
classes marital;
model dep_sx = marital / solution;
contrast 'married vs separated' marital -1 0 1 0 0 ;  * for question A;
contrast 'mean of widowed and divorced vs married' marital 2 -1 0 -1 0; * for question B;
run;


# 3.2, Linear regression
data temp;
set fhs.fhsdat;
proc sgplot data=temp;
scatter x=_age y=_bmi;
run;

proc glm data=temp2;
mode _bmi = _age agesq agecu / SS1;
title "cubic linear regression";
run;

proc reg data=temp3;
model _bmi = _age trr diastolc systolic / 
selection = stepwise sle=0.1 sls=0.1 details;
title "backward and forward selection";
plot residual. * predicted. ;
run;

# 3.2.2, mixed model (random effect and fixed effect)
proc mixed data=long;
class id group time;
model y = group time time*group / noint s ddfm=betwithin;
repeated time / subject=id type=un;

contrast 'trt-plc at month 1' group 1 -1 group*time 1 0 0 -1 0 0 /E;
contrast 'trt-plc at month 2' group 1 -1 group*time 0 1 0 0 -1 0 /E;
contrast 'trt-plc at month 3' group 1 -1 group*time 0 0 1 0 0 -1 /E;
run;


# 3.3, Logistic (GLM)
proc logistic data=hw order=internal;
model opentime = age sex nowdep haspets worknow voluntr / RL details 
selection=stepwise sle=.1 sls=.1 lackfit rsquare;
run;

data one;
agecoeff = 0.1255;	* by 1 increment;
sexcoeff = 0.4837;	*male for 1, female for 2;
workcoeff = -0.6096;	* work for 1, non work for 2;
voluncoeff = -0.3335;	* volun for 1, non volun for 0;
intercept = -11.2014;
logodds = intercept + 75*agecoeff + 1*sexcoeff + 0*workcoeff + 0*voluncoeff;
run;

# 3.3.2, Possion (GLM)
proc genmod data=melanoma order=data;
class age region;
model cases=region age / dist=poisson link=log offset=ltotal;
run;

# 3.4, Survival analysis 
data sum2;
length end_date 8;
event=0;	*event means censored;
format end_date mmddyy8.;
set sum;
if outdate=. & falldat=. then do; end_date='08Dec1990'd; event=1; end;
else if outdate=. & falldat>. then end_date=falldat;
else if outdate>. & falldat=. then do; end_date=outdate; event=1; end;
else if outdate < falldat then do; end_date=falldat; event=0; end;  *5 special situation;
else if outdate > falldat then end_date=falldat;
else end_date=.;
long = end_date - date;
if long < 0 then delete;
run;
proc contents data=sum2; *(1355 obs);
run;

proc lifetest method=life intervals=0 to 1200 by 50 plots=(s,h);
time long*event(1);
strata sex;
title 'homework 2';
run;

proc lifetest data=sa_11 plots=(s, ls, lls)
method=KM Nelson conftype=LOGLOG confband=all
plots=survival(cl cb=ep strata=panel);
time time*indicator(0);
strata group;
run;

# 4, study design
# 4.1, resampling
proc surveyselect data=hrt out=hrtboot_placebo
seed=12345678
method=urs
samprate=1
outhits
rep=100;
where (treatment="Placebo");
run;

proc surveyselect data=hrt out=hrtboot_aspirin
seed=12345678
method=urs
samprate=1
outhits
rep=100;
where (treatment="Aspirin");
run;


proc surveyselect data=jhs method=srs
seed=20170425 n=(15 12 13) out=sample2;
strata grade;
run;

# 4.2, study plan assignment
proc plan seed=&seed;
factors block=40 random group=6 / noprint;
output out=randschedule
group nvals=(1 1 1 2 2 2)
random;
run;

%macro strrand;
%let seed=20170328;

*generate 3 seeds save in dataset ranum4, age cat 1 2 3;

data rannum4(keep=seed);
do i=1 to 3;
rannum=ranuni(&seed);
seed=int(rannum*1000000);
end;
run;

proc sql noprint;
select seed 
into :seedsz1 - :seedsz3
from rannum4
quit;

%do i=1 %to 3;
proc plan seed=&&seedsz&i;
factors block=30 random group=6 /noprint;
output out=datasz&i
group nvals=(1 2 2 3 3 3)
random;
run;

data data&i;
set datasz&i;
age=&i;
run;
%end;

data sch;
set
%do i=1 %to 3;
data&i %end;
;
run;

proc print data=sch;
run;
%mend strrand;

%strrand;


# 4.3, power calculation
proc power;
twosamplesurvival test=logrank
hazardratio=0.576
refsurvexphazard=0.241
accrualtime=2
totaltime=5
ntotal=.
power=0.8;
run;

# 4.4, group sequential design
ods graphics on;
proc seqdesign altref=0.15
				errspend
				stopprob
				plots=errspend
				;

onesidePeto: design method=peto
			nstages=3
			alt=upper
			stop=both
			alpha=0.025 beta=0.1;
samplesize model=twosamplefreq (nullprop=0.6  test=prop);
ods output Boundary=Bnd_Count;
ods graphics off;





















