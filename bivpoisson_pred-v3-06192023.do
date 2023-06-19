*! version 1.0.1 , June-18-2023
*! Author: Abbie Zhang
*! Co-Authors: James Fisher, Joseph Terza
*! Website: https://github.com/zhangyl334/bivpoisson_predict
*! Support: 


*!***********************************************************!
*!     Count Valued Seemingly Unrelated Regression   	    *!
*!***********************************************************!


/* 1: DEEP PARAMETER ESTIMATION */


capture program drop bivpoisson
program define bivpoisson, sortpreserve eclass
	version 17.0
	/*Syntax: (t1 = dvars1) (t2 = dvars2) [if]
	Model automatically adds a constant. 
	Code is not-yet robust to syntax errors.*/
	
	/*Parsing of Inputs*/
	/*Equations and If Statement*/
		gettoken eq1 0 : 0, parse(") (") match(parns)
		gettoken eq2 0 : 0, parse(") (") match(parns)
		local ifstatement = "`0'"
	/*Dependent Variables from Equations*/
		gettoken dep1 eq1 : eq1, parse(" = ")
		gettoken dep2 eq2 : eq2, parse(" = ")
	/*Independent Variables from Equations*/		
		gettoken t eq1 : eq1, parse("=") /*Remove equals sign*/
		gettoken indep1 : eq1, parse("=")
		gettoken t eq2 : eq2, parse("=") /*Remove equals sign*/
		gettoken indep2 : eq2, parse("=")
	/*Parsed strings:
		dep1 = eq1 dependent variable
		dep2 = eq2 dependent variable
		indep1 = eq1 independent variables
		indep2 = eq2 independent variables
		ifstatement = sample to use for analysis
	*/
	
	/*Mark Sample*/
		tempvar touse
		qui gen byte `touse' = 0
		qui replace `touse' = 1 `ifstatement'	
		markout	`touse' `dep1' `dep2' `indep1' `indep2'	 /*Drop missing variables*/
	
	/*Check Variables*/
		/*Check that targets are positive, have 
		positive total variation, and are integer valued*/
		/*Eq 1*/ 
			qui sum `dep1' if `touse' == 1
			if r(min) < 0 {
				dis in green "{bf:`dep1'} is negative"
				exit 2000
			}	
			if r(Var) == 0 {
				dis in green "{bf:`dep1'} does not vary"
				exit 2000
			}	
			tempvar tmp
			qui gen `tmp' = (`dep1' - int(`dep1'))^2 if `touse' == 1
			qui sum `tmp' if `touse' == 1
			if r(sum) >0 {
				dis in green "{bf:`dep1'} is not integer valued"
				exit 2000
			}
		/*Eq2*/
			qui sum `dep2' if `touse' == 1
			if r(min) < 0 {
				dis in green "{bf:`dep2'} is negative"
				exit 2000
			}	
			if r(Var) == 0 {
				dis in green "{bf:`dep2'} does not vary"
				exit 2000
			}	
			tempvar tmp
			qui gen `tmp' = (`dep2' - int(`dep2'))^2 if `touse' == 1
			qui sum `tmp' if `touse' == 1
			if r(sum) > 0 {
				dis in green "{bf:`dep2'} is not integer valued"
				exit 2000
			}
		/*Check that target variables are not overly zero-inflated.
		Current test: poisson regression returns coefficient >= 0.
		Needs work*/
		/*Eq1*/
			qui: poisson `dep1' if `touse' == 1
			scalar tmp = e(b)[1,1]
			if tmp < 0 {
				dis in green "{bf:`dep1'} is zero-inflated"
				exit 2000
			}
		/*Eq2*/
			qui: poisson `dep2' if `touse' == 1
			scalar tmp = e(b)[1,1]
			if tmp < 0 {
				dis in green "{bf:`dep2'} is zero-inflated"
				exit 2000
			}
		/*Check for colinear feature variables and remove them*/
		/*Eq1*/
			qui _rmcoll `indep1' if `touse' == 1, forcedrop
			local indep1 "`r(varlist)'"
			if r(k_omitted) > 0 {
				dis in green "{bf:EQ1} several independent variables are colinear, automatically dropping them"
				dis in green "{bf:EQ1} revised independent variables are: `indep1'"
			}
		/*Eq2*/
			qui _rmcoll `indep2' if `touse' == 1, forcedrop
			local indep2 "`r(varlist)'"
			if r(k_omitted) > 0 {
				dis in green "{bf:EQ2} several independent variables are colinear, automatically dropping them"
				dis in green "{bf:EQ2} revised independent variables are: `indep2'"
			}
		
	/*Starting Values*/
	/*Eq1 via Poisson Regression*/
		qui poisson `dep1' `indep1' if `touse' == 1
		if _rc == 0{
			tempname cb1
			mat `cb1' = e(b)
			local ll_1 = e(ll)
		}
		if _rc !=0 {
			dis in green "{bf:EQ1} Initial values could not be estimated"
			exit 2000
		}
	/*Eq2 via Poisson Regression*/
		qui poisson `dep2' `indep2' if `touse' == 1
		if _rc == 0{
			tempname cb2
			mat `cb2' 	= e(b)
			local ll_2 = e(ll)
		}
		if _rc !=0 {
			dis in green "{bf:EQ2} Initial values could not be estimated"
			exit 2000
		}
	/*Starting Values for Rho via Assumption; needs work.*/
		tempname sigma12
		mat `sigma12' = (0)
		mat colnames `sigma12' = "/:sigma12"
	
	/*Mata Load*/
		/*Data
			We save into Y1, Y2, X1, and X2.  These will be overwritten. */
			tempvar cons
			qui gen `cons' = 1
			qui putmata Y1 = `dep1' if `touse' == 1, replace
			qui putmata Y2 = `dep2' if `touse' == 1, replace
			qui putmata X1 = (`indep1' `cons') if `touse' == 1, replace
			qui putmata X2 = (`indep2' `cons') if `touse' == 1, replace
		/*Initial Values*/
			mata: beta1_1n = st_matrix("`cb1'")
			mata: beta2_1n = st_matrix("`cb2'")
			mata: sigma12 = st_matrix("`sigma12'")
			mata: sigmasq1= 1  // Diagonal element of bi-NRV
			mata: sigmasq2 = 1  // See above.
		/*Parameters*/
			mata: quadpts = 30 // Number of quadrature points 
			mata: lims = (-5,5)  // Vector of numerical integration limits
			mata: limits = lims#J(rows(Y1),1,1) //Transformation
	
	/*Perform Estimation in Mata*/ 
		/*Setup Problem*/
			qui capture mata: mata drop BivPoissNorm
			qui mata: BivPoissNorm=moptimize_init()
			qui mata: moptimize_init_evaluator(BivPoissNorm, &BivPoissNormLF())
			qui mata: moptimize_init_evaluatortype(BivPoissNorm, "lf0")
			qui mata: moptimize_init_depvar(BivPoissNorm, 1, Y1)
			qui mata: moptimize_init_depvar(BivPoissNorm, 2, Y2)
			qui mata: moptimize_init_eq_indepvars(BivPoissNorm, 1, X1)
			qui mata: moptimize_init_eq_cons(BivPoissNorm, 1, "off") 
			qui mata: moptimize_init_eq_colnames(BivPoissNorm, 1, tokens("`indep1' _cons"))
			qui mata: moptimize_init_eq_indepvars(BivPoissNorm, 2, X2)
			qui mata: moptimize_init_eq_colnames(BivPoissNorm, 2, tokens("`indep2' _cons"))
			qui mata: moptimize_init_eq_cons(BivPoissNorm, 2,  "off" ) 
			qui mata: moptimize_init_eq_indepvars(BivPoissNorm, 3, "")
			qui mata: moptimize_init_eq_indepvars(BivPoissNorm, 4, "")
			qui mata: moptimize_init_eq_indepvars(BivPoissNorm, 5, "")
			qui mata: moptimize_init_eq_name(BivPoissNorm, 1, "`dep1'")
			qui mata: moptimize_init_eq_name(BivPoissNorm, 2, "`dep2'")
			qui mata: moptimize_init_eq_name(BivPoissNorm, 3, "sigmasq1")
			qui mata: moptimize_init_eq_name(BivPoissNorm, 4, "sigmasq2")
			qui mata: moptimize_init_eq_name(BivPoissNorm, 5, "sigma12")
		/*Initial Values*/
			qui mata: moptimize_init_eq_coefs(BivPoissNorm, 1, beta1_1n)
			qui mata: moptimize_init_eq_coefs(BivPoissNorm, 2, beta2_1n)
			qui mata: moptimize_init_eq_coefs(BivPoissNorm, 3, sigmasq1)
			qui mata: moptimize_init_eq_coefs(BivPoissNorm, 4, sigmasq2)
			qui mata: moptimize_init_eq_coefs(BivPoissNorm, 5, sigma12)
		/*Solve*/
			mata: moptimize(BivPoissNorm)
		/*Write results to console + Stata */
			mata: moptimize_result_display(BivPoissNorm)
			qui mata: moptimize_result_post(BivPoissNorm)
			/*Additional Entries for Ereturn*/
			ereturn local cmd "bivpoisson"
			ereturn local title "Bivariate Count Seemingly Unrelated Regression Estimation"
			ereturn local depvar1 `dep1'
			ereturn local indep1 `indep1'
			ereturn local depvar2 `dep2'
			ereturn local indep2 `indep2'
			ereturn local ifstatement "`ifstatement'"	
end

/*Mata Programs*/
	/*Quadrature Weights and Abscissa*/
	capture mata: mata drop GLQwtsandabs()
	mata 
	matrix GLQwtsandabs(real scalar quadpts)
	{
	  i = (1..quadpts-1)
	  b = i:/sqrt(4:*i:^2:-1) 
	  z1 = J(1,quadpts,0)
	  z2 = J(1,quadpts-1,0)
	  CM = ((z2',diag(b))\z1) + (z1\(diag(b),z2'))
	  V=.
	  ABS=.
	  symeigensystem(CM, V, ABS)
	  WTS = (2:* V':^2)[,1]
	  return(WTS,ABS') 
	} 
	end

	/*Integrand Bivariate Probit*/
	capture mata: mata drop BivPoissNormIntegrand()
	mata
	real matrix BivPoissNormIntegrand(real matrix xxu1, real matrix xxu2, /*
								   */ real matrix Y1, real matrix Y2, /*
								   */ real matrix xb1, real matrix xb2, /*
								   */ real matrix sigma12, real matrix sigmasq1, /*
								   */ real matrix sigmasq2)
	{
	lambda1=exp(xb1:+xxu1)
	lambda2=exp(xb2:+xxu2)
	
	poisspart=poissonp(lambda1,Y1):*poissonp(lambda2,Y2)
	
	SIGMA= sigmasq1,sigma12 \
           sigma12,sigmasq2
		   
	xxu=colshape(xxu1,1),colshape(xxu2,1)
	
	factor=rowsum((xxu*invsym(SIGMA)):*xxu)
	
	bivnormpart= (1:/(2:*pi():*sqrt(det(SIGMA))))/*
			   */ :*exp(-.5:*factor)
			   
	matbivnormpart=colshape(bivnormpart,cols(xxu1))
	
	integrandvals=poisspart:*matbivnormpart		 
	return(integrandvals)
	}	
	end

	/*2-D Integration Procedure*/
	capture mata: mata drop bivquadleg()
	mata
	real matrix bivquadleg(pointer(function) func, real matrix limits1, /*
						 */ real matrix limits2, real matrix wtsabs, /*
						 */ real matrix Y1, real matrix Y2, real matrix xb1, /*
						 */ real matrix xb2, real matrix sigma12, /*
						 */ real matrix sigmasq1, real matrix sigmasq2)
	{
	wts=wtsabs[.,1]'
	abcissae=wtsabs[.,2]'
	quadpts=rows(wtsabs)
	constant11=(limits1[.,2]:-limits1[.,1]):/2
	constant12=(limits1[.,2]:+limits1[.,1]):/2
	constant21=(limits2[.,2]:-limits2[.,1]):/2
	constant22=(limits2[.,2]:+limits2[.,1]):/2
	abcissaeC=J(1,quadpts,1)#abcissae'
	abcissaeR=abcissaeC'
	vecabcissaeC=rowshape(abcissaeC,1)
	vecabcissaeR=rowshape(abcissaeR,1)
	bigargs1=vecabcissaeC#constant11:+constant12
	bigargs2=vecabcissaeR#constant21:+constant22
	funvals=(*func)(bigargs1, bigargs2, Y1, Y2, xb1, xb2, sigma12, sigmasq1, sigmasq2)
	bigwts=wts'*wts
	vecbigwts=rowshape(bigwts,1)
	summand=constant11:*constant21:*(vecbigwts:*funvals)
	integapprox=colsum(summand')
	return(integapprox')
	}
	end

	/*Objective Function for Bivariate Probit*/
	capture mata: mata drop BivPoissNormLF()
	mata
	function BivPoissNormLF(transmorphic BivPoissNorm, real scalar todo, /*
						 */ real rowvector b, real matrix fv, real matrix SS, /*
						 */ real matrix HH) 
	{
	Y1 = moptimize_util_depvar(BivPoissNorm, 1)
	Y2 = moptimize_util_depvar(BivPoissNorm, 2)
	xb1 = moptimize_util_xb(BivPoissNorm, b, 1)
	xb2 = moptimize_util_xb(BivPoissNorm, b, 2)
	sigmasq1 = moptimize_util_xb(BivPoissNorm, b, 3)
	sigmasq2 = moptimize_util_xb(BivPoissNorm, b, 4)
	sigma12 = moptimize_util_xb(BivPoissNorm, b, 5)
	external quadpts
	external limits
	wtsandabs=GLQwtsandabs(quadpts)	
	likeval=bivquadleg(&BivPoissNormIntegrand(), limits, limits, wtsandabs,
						Y1, Y2, xb1, xb2, sigma12, sigmasq1, sigmasq2)     
	fv=ln(likeval)
	}
	end
	
	

	
	
*! version 1.1.0 , June-18-2023
*! Author: James Fisher, Abbie Zhang
*! Website: yileizhang.com
*! Support: zhangyl334@gmail.com


/***************************************************************************/
/*    			 bivpoisson prediction							           */
/* The goal for this command is to use the deep parameter estimates        */
/* produced by bivpoisson (beta1_hats, beta2_hats, and rho_hat             */
/* to predict the count valued outcomes Y1_hat and Y2_hat.                 */
/* Y1_hat = exp(X1*beta1_hat), Y2_hat = exp(X2*beta2_hat).                 */
/* The user can supply the value the policy variable (binary for now)      */
/* with the goal of predicting counterfactural value of Y1 (or Y2)         */

/* it is revised from rbiprobit, where there are 4 predicted probabilities */
/* 11, 10, 01 and 00. Revision, we predict based on conditional mean       */
/* function: exponentiated linear index llambda1 = exp(XoBo2 + XB1)        */
/* and llambda2 = exp(XoBo2 + XB2). */

/*Y1_1 (and Y2_1) represents the prediction of first (second) outcome when policy */
/* variable equal to 1,  Y1_0 (and Y2_0) represents the prediction of      */
/* first (second) outcome when policy variable equal to 0. (untreated)     */
/***************************************************************************/
/** Need to rewrite the functions and incorporate into the main program   **/
/** after parsing the ereturn and Xs. **/



/**
step 1) 
in Mata, parse the ereturn e(b) into beta1, beta2, and the 
rest ancilary parameters

estbeta1 = bbeta[1..3]  --> estbeta1 = `e(b)'[1..(`e(K)'-3)/2]
estbeta2 = bbeta[4..6]  --> estbeta2 = `e(b)'[(`e(K)'-3)/2+1..`e(K)'-3]
estsigma1 = bbeta[7]    --> estsigmasq1 = `e(b)'[`e(K)'-2]
estsigma2 = bbeta[8]    --> estsigmasq2 = `e(b)'[`e(K)'-1]
estsigma12 = bbeta[8]   --> estsigma12 = `e(b)'[`e(K)']

step 2)
parse the input X1, X2 vector into [X1o, X1, cons], and [X1o, X1, cons]
-- user input design: instruct users to put policy variable as the first variable 
User will write: 
Y1 (Xo=1,X1)
Y2 (Xo=1,X2)

or:
Y1 (Xo=0,X1)
Y2 (Xo=0,X2)

(xo is the policy var, same across equations, X1 and X2 can vary)

/*FURTHER Parsing of Inputs (how to seperate var1, and var2..... ?  */
	/*Equations and If Statement*/
		gettoken eq1 0 : 0, parse(") (") match(parns)
		gettoken eq2 0 : 0, parse(") (") match(parns)
		local ifstatement = "`0'"
	/*Dependent Variables from Equations*/
		gettoken dep1 eq1 : eq1, parse(" ( ")
		gettoken dep2 eq2 : eq2, parse(" ( ")
	/*Policy Variables from Equations*/		
		gettoken t eq1 : eq1, parse(" =") 
		gettoken policyvar : eq1, parse(" =")
		gettoken t eq2 : eq2, parse(" =") 
		gettoken policyvar : eq2, parse(" =")
	/*Independent Variables from Equations*/		
		gettoken t eq1 : eq1, parse(",") 
		gettoken indep1 : eq1, parse("=")
		gettoken t eq2 : eq2, parse("=,") 
		gettoken indep2 : eq2, parse("=")
		
step 3: calculate the exp(linear index): potential outcomes

qui mata etamat = bivariate normal 

predicted outcomes (when user input provides the policy variables value = 1 or 0)

qui mata: LAMBDA_y1=exp(tokens("`policyvar')*estbeta1[.,1]' :+ tokens("`indep1' _cons")*estbeta1[.,2..col(estbeta1)]':+sqrt(estsigmasq1):*etamat[.,1]) 

qui mata: LAMBDA_y2=exp(tokens("`policyvar')*estbeta2[.,1]' :+ tokens("`indep2' _cons")*estbeta2[.,2..col(estbeta2)]':+sqrt(estsigmasq2):*etamat[.,1]) 


POPULATION ATE (potential outcomes of the population when the treatment is 
changed from untreated to treated:
ATE (Y1):

qui mata: LAMBDA1_DELTA=exp(tokens("`1')*estbeta1[.,1]' :+ tokens("`indep1' _cons")*estbeta1[.,2..col(estbeta1)]':+sqrt(estsigmasq1):*etamat[.,1]) 

qui mata: LAMBDA1=exp(tokens("`0')*estbeta1[.,1]' :+ tokens("`indep1' _cons")*estbeta1[.,2..col(estbeta1)]':+sqrt(estsigmasq1):*etamat[.,1]) 

ATE (Y2):
qui mata: LAMBDA2_DELTA=exp(tokens("`1')*estbeta2[.,1]' :+ tokens("`indep2' _cons")*estbeta2[.,2..col(estbeta2)]':+sqrt(estsigmasq2):*etamat[.,1]) 

qui mata: LAMBDA2=exp(tokens("`0')*estbeta2[.,1]' :+ tokens("`indep2' _cons")*estbeta2[.,2..col(estbeta2)]':+sqrt(estsigmasq2):*etamat[.,1]) 


qui mata: ate_y1 = MEANCOUNT(LAMBDA1_DELTA):- MEANCOUNT(LAMBDA1)

qui mata: ate_y2 = MEANCOUNT(LAMBDA2_DELTA):- MEANCOUNT(LAMBDA2)

qui mata: AATE_y1=mean(ate_y1)
qui mata: AATE_y2=mean(ate_y2)


AAIE=mean(aie)

/*SAVE all Results for output */
			ereturn local AATE_y1 `AATE_y1'
			ereturn local Y1_predict `LAMBDA_y1'
			ereturn local AATE_y2 `AATE_y2'
			ereturn local Y2_predict `LAMBDA_y1'
	
provide a short description at the end of the estimate:
(find syntax for it)
Y1_predict report the population potential outcome of Y1
when the policy varibale is set to be `xo'. 
Y2_predict report the population potential outcome of Y2
when the policy varibale is set to be `xo'. 
**/


/*
keep the mean count function. 
Do not use the AIEfunction, but incorporate it into the final ATE estimate.
**/

****************************************************************
/***************************************************************
    bivpoisson_predict Sytax: (Y1 (Xo=1,X1)) (Y2 (Xo=1,X2)) [if]
	Or keep it the same as bivpoisson?
	Model automatically adds a constant. 	
	Defaul Output: 
	1) potential outcomes Y1, Y2 (when whole population is treated and 
	   untreated respectively)
	2) ATE of Y1, Y2 (Population Average Treatment Effects when 
	  policy variable of interest: xo is alterded from 0 to 1. )
	3) conditional potential outcomes: Y1 | (xo=user defined value)
	                               and Y2 | (xo=user defined value)

	recall: bivpoisson Syntax: (t1 = dvars1) (t2 = dvars2) [if]
***************************************************************/

capture program drop bivpoisson_predict
program define bivpoisson_predict, sortpreserve eclass
	version 17.0
	
/* Parsing of Inputs (how to seperate var1, and var2..... ?  */
	/*Equations and If Statement*/
		gettoken eq1 0 : 0, parse(") (") match(parns)
		gettoken eq2 0 : 0, parse(") (") match(parns)
		local ifstatement = "`0'"
	/*Dependent Variables from Equations*/
		gettoken dep1 eq1 : eq1, parse(" ( ")
		gettoken dep2 eq2 : eq2, parse(" ( ")
	/*Policy Variables from Equations*/		
		gettoken t eq1 : eq1, parse(" =") 
		gettoken policyvar : eq1, parse(" =")
		gettoken t eq2 : eq2, parse(" =") 
		gettoken policyvar : eq2, parse(" =")
	/*Independent Variables from Equations*/		
		gettoken t eq1 : eq1, parse(",") 
		gettoken indep1 : eq1, parse("=")
		gettoken t eq2 : eq2, parse("=,") 
		gettoken indep2 : eq2, parse("=")	
	
/*************************************************************
Call the bivpoisson command to get deep parameter estimates
***************************************************************
user should install bivpoisson by typing: ssc install bivpoisson 
in Stata demand line (Stata 17 or higher)
**************************************************************/
ssc install bivpoisson

bivpoisson (`dep1' = tokens("`policyvar' `indep1' _cons")) (`dep2' = tokens("`policyvar' `indep2' _cons"))

ereturn local title "Bivariate Count Seemingly Unrelated Regression Estimation"
			ereturn local depvar1 `dep1'
			ereturn local indep1 `indep1'
			ereturn local depvar2 `dep2'
			ereturn local indep2 `indep2'
			ereturn local ifstatement "`ifstatement'"	
			
qui mata estbeta1 = `e(b)'[1..(`e(K)'-3)/2]
qui mata estbeta2 = `e(b)'[(`e(K)'-3)/2+1..`e(K)'-3]
qui mata estsigmasq1 = `e(b)'[`e(K)'-2]
qui mata estsigmasq2 = `e(b)'[`e(K)'-1]
qui mata estsigma12 = `e(b)'[`e(K)']

/** calculate the exp(linear index): potential outcomes **/
/** Condtional Mean Function Y1 = exp(Xo*betao + X1*estbeta1 + sigma1*eta1). */
/** Similar for Y2. **/

/**************************************************************************
User will need to supply policy variable value of interest Xo=1 or 0, 
                          and covariates: X1, X2. **
**************************************************************************/
qui mata: SIGMA = 1, `estsigma12' \
                 `estsigma12', 1			 
qui mata: transmat=cholesky(SIGMA)
qui mata etamat=(transmat*(rnormal(1, rows(`Y1'),0,1) \ rnormal(1,rows(`Y1'),0,1)))'

qui mata: LAMBDA_y1=exp(tokens("`policyvar')*estbeta1[.,1]' :+ tokens("`indep1' _cons")*estbeta1[.,2..col(estbeta1)]':+sqrt(estsigmasq1):*etamat[.,1]) 

qui mata: LAMBDA_y2=exp(tokens("`policyvar')*estbeta2[.,1]' :+ tokens("`indep2' _cons")*estbeta2[.,2..col(estbeta2)]':+sqrt(estsigmasq2):*etamat[.,1]) 


/***************************************************************************
POPULATION ATE -- ATE (Y1):
 (potential outcomes of the population when the treatment is 
conterfacturally changed from untreated to treated: 
***************************************************************************/

qui mata: LAMBDA1_DELTA=exp(tokens("`1')*estbeta1[.,1]' :+ tokens("`indep1' _cons")*estbeta1[.,2..col(estbeta1)]':+sqrt(estsigmasq1):*etamat[.,1]) 

qui mata: LAMBDA1=exp(tokens("`0')*estbeta1[.,1]' :+ tokens("`indep1' _cons")*estbeta1[.,2..col(estbeta1)]':+sqrt(estsigmasq1):*etamat[.,1]) 

/**
ATE (Y2):
**/

qui mata: LAMBDA2_DELTA=exp(tokens("`1')*estbeta2[.,1]' :+ tokens("`indep2' _cons")*estbeta2[.,2..col(estbeta2)]':+sqrt(estsigmasq2):*etamat[.,1]) 

qui mata: LAMBDA2=exp(tokens("`0')*estbeta2[.,1]' :+ tokens("`indep2' _cons")*estbeta2[.,2..col(estbeta2)]':+sqrt(estsigmasq2):*etamat[.,1]) 

qui mata: function MEANCOUNT(llambda){
COMmean=llambda
return(COMmean)
}

qui mata: ate_y1 = MEANCOUNT(LAMBDA1_DELTA):- MEANCOUNT(LAMBDA1)
qui mata: ate_y2 = MEANCOUNT(LAMBDA2_DELTA):- MEANCOUNT(LAMBDA2)
qui mata: AATE_y1=mean(ate_y1)
qui mata: AATE_y2=mean(ate_y2)
qui mata: AAIE=mean(aie)

/*SAVE all Results for output */ 
        ereturn local LAMBDA1_DELTA `predict_y1_all_treated'
		ereturn local LAMBDA1 `predict_y1_all_untreated'
		ereturn local LAMBDA2_DELTA `predict_y2_all_treated'
		ereturn local LAMBDA2 `predict_y2_all_untreated'
		ereturn local AATE_y1 `AATE_y1'
		ereturn local AATE_y2 `AATE_y2'
		ereturn local Y1_predict `predict_y1(PolicyVar=`policyvar')'
		ereturn local Y2_predict `predict_y2(PolicyVar=`policyvar')'
	

	   local LAMBDA1_DELTA "Prediction of `dep1' when all treated"
	   local LAMBDA1 "Prediction of `dep1 when all untreated"
	   	
	   local LAMBDA2_DELTA "Prediction of `dep2' when all treated"
	   local LAMBDA2 "Prediction of `dep2 when all untreated"
	  
	   local AATE_y1 "Average Treatment Effects of `PolicyVar' on `dep1'"
	   local AATE_y2 "Average Treatment Effects of `PolicyVar' on `dep2'"
	   
	   local Y1_predict "Prediction of `dep2' when PolicyVar=`policyvar'"
	   local Y2_predict "Prediction of `dep2' when PolicyVar=`policyvar'"

	   
	   end
**/
/**
STOPPD HERE. 06/19/2023 2:21 pm
***/

	

mata 

function MEANCOUNT(llambda,JJstar){

COMmean=llambda
return(COMmean)
}



real matrix AIEfunct(bbeta, AAIE)
{

external aie
external XD
external X
external etamat

estbeta1 = bbeta[1..3]
estbeta2 = bbeta[4..6]
estsigma1 = bbeta[7]
estsigma2 = bbeta[8]


XDBETA=XD*estbeta1[.,1..2]':+sqrt(estsigma1):*etamat[.,1]
XB=X*estbeta1[.,1..2]':+sqrt(estsigma1):*etamat[.,1]

LAMBDADELTA=exp(XDBETA)
LAMBDA=exp(XB)


aie = MEANCOUNT(LAMBDADELTA,Jstar):- MEANCOUNT(LAMBDA,Jstar)
AAIE=mean(aie)

return(AAIE)

}

end 

/****************************************************************/


program define bivpoisson_p, eclass
	
	version 11
	syntax [anything] [if] [in] [, SCores *]
	
	if ("`e(cmd)'" != "bivpoisson"){
		error 301
		dis in red "bivpoisson was not the last command"
	}
	
	if `"`scores'"' != ""{
		ml score `0'	
		exit
	}
	
	local myopts Y1_1 Y2_1 Y1_0 Y2_0

	_pred_se "`myopts'" `0'
	
	if (`s(done)') exit			
	local vtyp 	`s(typ)'			
	local varn 	`s(varn)'		
	local 0	`"`s(rest)'"'			
	

	*!	parse predict	
	syntax [if] [in] [, `myopts' noOFFset]
	
	local type `Y1_1'`Y2_1'`Y1_0'`Y2_0'
	
	tokenize `e(depvar)'
	local dep1 `1'
	local dep2 `2'
	

	tsunab dep1: `dep1'
	tsunab dep2: `dep2'
				
	rmTS `dep1'
	confirm variable `r(rmTS)'
	local dep1n 	`r(rmTS)'
	
	rmTS `dep2'
	confirm variable `r(rmTS)'
	local dep2n 	`r(rmTS)'
	
	*! parse the deep paramter estimates into 5 pieces:
	*! beta1_hat, beta2_hat, and 3 ancilary parameters
	tokenize `e(b)'
	local beta1_hat `1..(e(rank)-3)/2'
	local beta2_hat `(e(rank)-3)/2+1..(e(rank)-3)'
	local sigmasq1 `(e(rank)-2)'
	local sigmasq2 `(e(rank)-1)'
	local sigma12 `e(rank)'
	
	
	mata
	DELTA = 1
	DELTA_minus = -1
	
	external XD_plus
	external XD_minus

	external X
	external etamat
	external DELTA
	external DELTA_minus
	
	XD_plus=X
    XD_plus[.,1]=X[.,1]:+DELTA
	
	XD_minus = X
	XD_minus[.,1]=X[.,1]:+DELTA_minus

	
	/* Here we formulate the first type of prediction: everyone is treated*/
	/* Xo = 1 for all observations, we calculate the values of Y1_1, and Y2_1 */
	
	/* XB1_plus1 is the first-element-altered linear index of outcome 1  */
	/* XB1 is the orignal linear index of outcome 1 (using beta_hat).    */
	
	local XB1_plus1 =XD_plus*`beta1_hat'':+sqrt(`sigmasq1'):*etamat[.,1]
    local XB1       =X *`beta1_hat'':+sqrt(`sigmasq1'):*etamat[.,1]

	local XB2_plus1 =XD_plus*`beta1_hat'':+sqrt(`sigmasq2'):*etamat[.,1]
	local XB2       =X*`beta2_hat'':+sqrt(`sigmasq2'):*etamat[.,1]

	
	/* XB2_plus1 is the first-element-altered linear index of outcome 2  */
	/* XB2 is the orignal linear index of outcome 2 (using beta_hat).    */
	local XB1_plus1 =XD_plus*`beta2_hat'':+sqrt(`sigmasq1'):*etamat[.,1]
    local XB1       =X *`beta1_hat'':+sqrt(`sigmasq1'):*etamat[.,1]

	local XB2_plus1 =XD_plus*`beta2_hat'':+sqrt(`sigmasq2'):*etamat[.,1]
	local XB2       =X*`beta2_hat'':+sqrt(`sigmasq2'):*etamat[.,1]
	
	
	/* Now, we formulate the second type of prediction: everyone is un-treated*/
	/* Xo = 0 for all observations, we calculate the values of Y1_0, and Y2_0  */
	
	/* XB1_minus1 is the first-element-altered (by -1) linear index of outcome 1  */
	/* XB2_minus1 is the first-element-altered (by -1) linear index of outcome 2  */
	
	local XB1_minus1 =XD_minus*`beta1_hat'':+sqrt(`sigmasq1'):*etamat[.,1]

	local XB2_minus1 =XD_minus*`beta2_hat'':+sqrt(`sigmasq2'):*etamat[.,1]

	Y1_1 = exp(`XB1_plus1')
	Y2_1 = exp(`XB2_plus1')

	Y1_0 = exp(`XB1_minus1')
	Y2_0 = exp(`XB2_minus1')

    end

	
	    *!	predict first outcome conditional on policy var = 1.
	if "`type'" == "Y1_1"{
		local Y1_1 "Prediction of `dep1' when all treated"
		
		exit
	}	
		
		*!	predict second outcome conditional on policy var = 1.
	if "`type'" == "Y2_1"{
		local Y2_1 "Prediction of `dep2 when all treated"
		
		exit
	}	
		
	
	     *!	predict first outcome conditional on policy var = 0.
	if "`type'" == "Y2_0"{
		local Y1_0 "Prediction of `dep1 when all un-treated"
		
		exit
	}	
		
		
		 *!	predict second outcome conditional on policy var = 0.
	if "`type'" == "Y2_0"{
		local Y2_0 "Prediction of `dep2 when all un-treated"
		
		exit
	}	
		

	}	
	

	error 198
end



program define rmTS, rclass
	
	local tsnm = cond( match("`0'", "*.*"),  		/*
			*/ bsubstr("`0'", 			/*
			*/	  (index("`0'",".")+1),.),     	/*
			*/ "`0'")

	return local rmTS `tsnm'
end
	
	