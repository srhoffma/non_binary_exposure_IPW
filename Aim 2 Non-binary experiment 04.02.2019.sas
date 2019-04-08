/*********************************************************************************
*
*	Study of Environment Lifestyle & Fibroids (SELF)
*	Hoffman Dissertation -- Aim 2 Non-binary IPW Experiment
*	03/13/2019
*
*********************************************************************************/

title ' '; title2 ' '; title3 ' '; title4 " "; title5 " "; title6 " "; 

/*********************************************************************************
*	Non-binary Exposure IPW (Changed 3rd row confounder to 0)
*********************************************************************************/

data IPW_1;
input ID Exposure Outcome Confounder;
CARDS;
1	0	1	0
3	0	0	0
2	0	1	0
4	0	1	1
17	0	1	1
5	1	1	0
18	1	0	0
7	1	0	1
8	1	0	1
6	1	1	1
10	2	1	0
11	2	1	0
12	2	0	0
19	2	1	0
9	2	0	1
20	3	1	0
13	3	1	1
14	3	1	1
15	3	0	1
16	3	1	1
;
run;

ods rtf style=htmlblue;

proc freq data=IPW_1; tables Exposure Outcome Confounder; run;

	/*Step 1: Crude Model*/
	title3 "Crude (Unadjusted) Model";
	proc genmod data=IPW_1 desc;
	class Exposure (order = internal ref=first);
	model Outcome=Exposure/link=log dist=binomial;
	estimate "Risk of outcome by exposure 1 v 0" Exposure 1 0 0 -1;
	estimate "Risk of outcome by exposure 2 v 0" Exposure 0 1 0 -1;
	estimate "Risk of outcome by exposure 3 v 0" Exposure 0 0 1 -1;
	run;

	/*Step 2: PS Model*/
	title3 "Propensity Score Model";
	proc logistic data=IPW_1 descending;
	class Exposure (order = internal ref=last) Confounder;
	model Exposure = Confounder /link=glogit;
	output out=ps_data_1 predicted=ps;
	run;

	/*Step 3: PS Curves*/
	data ps_data_2;
	set ps_data_1;
	if Exposure = 0 then ps_0 = ps;
		else ps_0 = .; 
	if Exposure = 1 then ps_1 = ps;
		else ps_1 = .; 
	if Exposure = 2 then ps_2 = ps;
		else ps_2 = .; 
	if Exposure = 3 then ps_3 = ps;
		else ps_3 = .; 
	run;

	proc freq data=ps_data_2; tables ps_0 ps_1; run;

	%macro level (trt = );
	proc kde data=ps_data_2;
		univar ps_0 ps_1 ps_2 ps_3 / plots=densityoverlay;
		title4 "Propensity score distributions by treatment group";
		title5 "Multinomial Logistic Model PS";
		title6 "Probability of being in treatment group &trt";
		where _level_ = &trt;
	run;
	%mend;

	%level (trt = 0);
	%level (trt = 1);
	%level (trt = 2);
	%level (trt = 3);

	title4 " "; title5 " "; title6 " "; 

	/*Step 4: IPW*/

	/*Estimation of numerator of stabalized ip weights*/
	title3 "Empty Model For Stabilized Weight Numerators";
	proc logistic data=IPW_1 descending;
		ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
		class Exposure (order = internal ref=first);
		model Exposure = /link=glogit;
		output out=ps_data_1_n (keep= id _level_ pn_Exposure) p=pn_Exposure;
	run;

	/*Transposing the data so that numerators and denominators for weights are in one place*/
	data ps_data_1x; set ps_data_1; keep id Exposure _level_ ps; run;
	data for_set; set IPW_1; run;
	proc sort data=ps_data_1x; by id; run;
	proc transpose data=ps_data_1x out = ps_data_1y prefix = ps_; var ps; id _level_; by id; run;	
	proc sort data=ps_data_1_n; by id; run;
	proc transpose data=ps_data_1_n out = ps_data_1n_y prefix = mp_; var pn_Exposure; id _level_; by id; run;
	proc sort data=for_set; by id; run;
	proc sort data=ps_data_1y; by id; run;
	proc sort data=ps_data_1n_y; by id; run;

	data ps_data_3;
	merge for_set ps_data_1y ps_data_1n_y;
	by id;
	/*Inverse probability of actual exposure*/
	if Exposure=0 then ipw= 1/ps_0;
	if Exposure=1 then ipw= 1/ps_1;
	if Exposure=2 then ipw= 1/ps_2;
	if Exposure=3 then ipw= 1/ps_3;
	/*Stabilized ipw: replacing numerator with marginal probability of actual exposure*/
	if Exposure=0 then s_ipw= mp_0/ps_0;
	if Exposure=1 then s_ipw= mp_1/ps_1;
	if Exposure=2 then s_ipw= mp_2/ps_2;
	if Exposure=3 then s_ipw= mp_3/ps_3;
	run;

	title3 "Distributions Of Propensity Scores and Weights";
	proc means data=ps_data_3 mean; 
	class Exposure; 
	var ipw s_ipw;
	run;

	proc univariate data=ps_data_3;
	var ipw s_ipw;
	id id;
	run;

	/*Step 5: Table 1*/ /*LEFT OFF HERE*/ /*LEFT OFF HERE*/ /*LEFT OFF HERE*/ /*LEFT OFF HERE*/

	/*Table 1 - Unweighted*/
	title3 "Unweighted Covariate Balance Table";
	proc means data=ps_data_3 noprint; 
	class Exposure; 
	var Confounder;
	output out = new (drop=_type_ _freq_) mean = ;
	run;

	proc transpose data=new out = new2 (drop=_label_) prefix = pct_; id Exposure; run;
	proc sort data=new2; by _name_; run;

	proc means data=ps_data_3 noprint; 
	class Exposure; 
	var Confounder;
	output out = new3 (drop=_type_ _freq_) std = ;
	run;

	proc transpose data=new3 out = new4 (drop=_label_) prefix = std_; id Exposure; run;
	proc sort data=new4; by _name_; run;

	data new5;
	merge new2 new4;
	by _name_;
	mean_abs_std_diff = round(mean(
		abs((pct_1 - pct_0)/sqrt(((std_1*std_1) + (std_0*std_0))/2)),
		abs((pct_1 - pct_2)/sqrt(((std_1*std_1) + (std_2*std_2))/2)),
		abs((pct_1 - pct_3)/sqrt(((std_1*std_1) + (std_3*std_3))/2)),
		abs((pct_2 - pct_0)/sqrt(((std_2*std_2) + (std_0*std_0))/2)),
		abs((pct_2 - pct_3)/sqrt(((std_2*std_2) + (std_3*std_3))/2)),
		abs((pct_3 - pct_0)/sqrt(((std_3*std_3) + (std_0*std_0))/2))	
			), .001);
	drop std_0-std_3;
	run;

	proc print data=new5; run;

	/*Table 1 - IPW Weighted*/
	title3 "IPW Weighted Covariate Balance Table";
	proc means data=ps_data_3 noprint; 
	class Exposure; 
	var Confounder;
	weight ipw;
	output out = new (drop=_type_ _freq_) mean = ;
	run;

	proc transpose data=new out = new2 (drop=_label_) prefix = pct_; id Exposure; run;
	proc sort data=new2; by _name_; run;

	proc means data=ps_data_3 noprint; 
	class Exposure; 
	var Confounder;
	weight ipw;
	output out = new3 (drop=_type_ _freq_) std = ;
	run;

	proc transpose data=new3 out = new4 (drop=_label_) prefix = std_; id Exposure; run;
	proc sort data=new4; by _name_; run;

	data new5;
	merge new2 new4;
	by _name_;
	mean_abs_std_diff = round(mean(
		abs((pct_1 - pct_0)/sqrt(((std_1*std_1) + (std_0*std_0))/2)),
		abs((pct_1 - pct_2)/sqrt(((std_1*std_1) + (std_2*std_2))/2)),
		abs((pct_1 - pct_3)/sqrt(((std_1*std_1) + (std_3*std_3))/2)),
		abs((pct_2 - pct_0)/sqrt(((std_2*std_2) + (std_0*std_0))/2)),
		abs((pct_2 - pct_3)/sqrt(((std_2*std_2) + (std_3*std_3))/2)),
		abs((pct_3 - pct_0)/sqrt(((std_3*std_3) + (std_0*std_0))/2))	
			), .001);
	drop std_0-std_3;
	run;

	proc print data=new5; run;

	/*Table 1 - Stabilized IPW Weighted*/
	title3 "Stabilized IPW Weighted Covariate Balance Table";
	proc means data=ps_data_3 noprint; 
	class Exposure; 
	var Confounder;
	weight s_ipw;
	output out = new (drop=_type_ _freq_) mean = ;
	run;

	proc transpose data=new out = new2 (drop=_label_) prefix = pct_; id Exposure; run;
	proc sort data=new2; by _name_; run;

	proc means data=ps_data_3 noprint; 
	class Exposure; 
	var Confounder;
	weight s_ipw;
	output out = new3 (drop=_type_ _freq_) std = ;
	run;

	proc transpose data=new3 out = new4 (drop=_label_) prefix = std_; id Exposure; run;
	proc sort data=new4; by _name_; run;

	data new5;
	merge new2 new4;
	by _name_;
	mean_abs_std_diff = round(mean(
		abs((pct_1 - pct_0)/sqrt(((std_1*std_1) + (std_0*std_0))/2)),
		abs((pct_1 - pct_2)/sqrt(((std_1*std_1) + (std_2*std_2))/2)),
		abs((pct_1 - pct_3)/sqrt(((std_1*std_1) + (std_3*std_3))/2)),
		abs((pct_2 - pct_0)/sqrt(((std_2*std_2) + (std_0*std_0))/2)),
		abs((pct_2 - pct_3)/sqrt(((std_2*std_2) + (std_3*std_3))/2)),
		abs((pct_3 - pct_0)/sqrt(((std_3*std_3) + (std_0*std_0))/2))	
			), .001);
	drop std_0-std_3;
	run;

	proc print data=new5; run;

	title3 "";

	/*Step 6: Final Model*/
	title3 "IPW Weighted Model";
	proc genmod data=ps_data_3 desc;
	class id Exposure (order = internal ref=first);
	weight ipw;
	model Outcome=Exposure/link=log dist=binomial;
	repeated subject = id /type = ind;
	estimate "Risk of outcome by exposure 1 v 0" Exposure 1 0 0 -1;
	estimate "Risk of outcome by exposure 2 v 0" Exposure 0 1 0 -1;
	estimate "Risk of outcome by exposure 3 v 0" Exposure 0 0 1 -1;
	run;

	proc freq data=ps_data_3;
	weight ipw;
	table Outcome*Exposure /chisq;
	run;

	title3 "Stabilized IPW Weighted Model";
	proc genmod data=ps_data_3 desc;
	class id Exposure (order = internal ref=first);
	weight s_ipw;
	model Outcome=Exposure/link=log dist=binomial;
	repeated subject = id /type = ind;
	estimate "Risk of outcome by exposure 1 v 0" Exposure 1 0 0 -1;
	estimate "Risk of outcome by exposure 2 v 0" Exposure 0 1 0 -1;
	estimate "Risk of outcome by exposure 3 v 0" Exposure 0 0 1 -1;
	run;

	/*Step 7: Traditional Model for Comparison*/
	title3 "Traditional Log Binomial Multivariate Regression";
	proc genmod data=IPW_1 desc;
	class Exposure (order = internal ref=first);
	model Outcome=Exposure Confounder/link=log dist=binomial;
	estimate "Risk of outcome by exposure 1 v 0" Exposure 1 0 0 -1;
	estimate "Risk of outcome by exposure 2 v 0" Exposure 0 1 0 -1;
	estimate "Risk of outcome by exposure 3 v 0" Exposure 0 0 1 -1;
	run;

	title ' '; title2 ' '; title3 ' '; title4 " "; title5 " "; title6 " "; 

	ods rtf close; 

/*********************************************************************************
*	Diagnostics
*********************************************************************************/

	/*Note: S_IPW and IPW weights matched what I obtained by hand a priori in Excel*/

	/*Sum of Weights, 20 and 80 respectively, 1x and 4x respectively*/
	/*S_IPW weights add up to the sample size, while non-standardized IPW 
		weights add up to t * sample size, where t = number of trt groups*/
	/*proc sql;
	select exposure, sum(s_ipw) as s_ipw_sum, sum(ipw) as ipw_sum
	from ps_data_3
	group by exposure;
	quit;

	/*Still, within each treatment group, the confounder distribution looks like
		that of the total sample*/
	/*Does the confounder distribution within each exposure group look like the 
		distribution in the total sample? Yes!*/
	/*proc means data=ps_data_3; var confounder; title 'Unweighted'; run;
	proc means data=ps_data_3; class exposure; var confounder; title 'Unweighted'; run;
	proc means data=ps_data_3; class exposure; var confounder; weight s_ipw; title 'S_IPW weighted'; run;

	
	proc freq data=ps_data_3; tables confounder; title 'Unweighted'; run;
	proc freq data=ps_data_3; tables exposure*confounder; title 'Unweighted'; run;
	proc freq data=ps_data_3; tables exposure*confounder; weight s_ipw; title 'S_IPW weighted'; run;

	/*Does this still hold true when just looking at two groups, such as
		comparing trt = 3 to trt = 0?*/
	/*proc means data=ps_data_3; var confounder; where exposure in (0,3); title 'Unweighted'; run;
	proc means data=ps_data_3; var confounder; weight s_ipw; where exposure in (0,3); title 'S_IPW weighted'; run;

	ods rtf close;

	proc datasets library=work kill;run; quit;

/*********************************************************************************
*	End of Code
*********************************************************************************/
