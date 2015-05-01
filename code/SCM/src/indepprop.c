
#include <math.h>
#include <unistd.h>

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#include "SCM.h"

typedef struct R_SCM_matrix_t {
  double *M;
  int nrow, ncol;
} R_SCM_matrix_t;

typedef struct R_SCM_imatrix_t {
  int *M;
  int nrow, ncol;
} R_SCM_imatrix_t;

#define MATRIX(A, i, j) (((A).M)[(A).nrow*(j)+(i)])

SEXP R_SCM_getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

double R_SCM_getNumericElement(SEXP vector, const char *str) {
  SEXP names=getAttrib(vector, R_NamesSymbol);
  double elmt=0.0;
  int i;

  for (i = 0; i < length(vector); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = REAL(vector)[i];
      break;
    }
  return elmt;
}

double R_SCM_target_dens(double x, double eta0, double eta1, double mu, 
			 double sig2, double thresh) {
  if (x < thresh) {
    return dnorm(x, mu, sqrt(sig2), /*give_log=*/ 0);
  } else {
    return 1.0 / (1.0 + exp(-(eta0+eta1*x))) * dnorm(x, mu, sqrt(sig2), 
						 /*give_log=*/ 0);
  }
}

double R_SCM_target_hessian(double x, double eta0, double eta1, 
			    double sig2) {
  double eta12=eta1*eta1;
  double tmp=1+exp(eta0+eta1*x);
  return -eta12 / tmp + eta12/(tmp*tmp) - 1.0/sig2;
}

typedef struct R_SCM_optim_data_t {
  double e0, e1, mu, sig2, thr;
} R_SCM_optim_data_t;

double R_SCM_fmin_callback(double x, void *extra) {
  R_SCM_optim_data_t *data=(R_SCM_optim_data_t*) extra;
  return -R_SCM_target_dens(x, data->e0, data->e1, data->mu, data->sig2,
			    data->thr);
}

/* 
 * The proposal distribution is a Normal, with its mean set as
 * the mode of the distribution we want to sample from
 * The variance should be set from the Fisher information 
 */

SEXP R_SCM_drawMissing(SEXP pstate, SEXP pmissIdx, SEXP peta0, SEXP peta1,
		       SEXP pthresholds, SEXP pdataMean, SEXP pdataVar,
		       SEXP ptopresR, SEXP ptopresE, SEXP ptopresT, SEXP pnoAcc,
		       SEXP pnoRej, SEXP pkj, SEXP ptj) {

	SEXP result, names;
	SEXP pX, pR, pE, pT;

	SEXP oX=R_SCM_getListElement(pstate, "X");
	SEXP oR=R_SCM_getListElement(pstate, "R");
	SEXP oE=R_SCM_getListElement(pstate, "E");
	SEXP oT=R_SCM_getListElement(pstate, "T");

	int Xrow=INTEGER(GET_DIM(oX))[0];
	int Xcol=INTEGER(GET_DIM(oX))[1];
	int Rrow=INTEGER(GET_DIM(oR))[0];
	int Rcol=INTEGER(GET_DIM(oR))[1];
	int Erow=INTEGER(GET_DIM(oE))[0];
	int Ecol=INTEGER(GET_DIM(oE))[1];
	int Trow=INTEGER(GET_DIM(oT))[0];
	int Tcol=INTEGER(GET_DIM(oT))[1];

	R_SCM_matrix_t X={ 0, Xrow, Xcol };
	R_SCM_matrix_t R={ 0, Rrow, Rcol };
	R_SCM_matrix_t E={ 0, Erow, Ecol };
	R_SCM_matrix_t T={ 0, Trow, Tcol };

	double *eta0=REAL(peta0);
	double *eta1=REAL(peta1);
	double *thresholds=REAL(pthresholds);
	double *mdataMean=REAL(pdataMean);
	double *dataVar=REAL(pdataVar);
	R_SCM_matrix_t dataMean={ mdataMean, Xrow, Xcol };

	int i, nomiss=INTEGER(GET_DIM(pmissIdx))[0];
	R_SCM_imatrix_t missIdx={ INTEGER(pmissIdx), nomiss, 2 };

	int noAcc=INTEGER(pnoAcc)[0], noRej=INTEGER(pnoRej)[0];

	int *topresR=LOGICAL(ptopresR);
	int *topresE=LOGICAL(ptopresE);
	int *topresT=LOGICAL(ptopresT);
	int *kj=INTEGER(pkj);
	int *tj=INTEGER(ptj);

  GetRNGstate();
	
	PROTECT(pX=duplicate(R_SCM_getListElement(pstate, "X")));
	PROTECT(pR=duplicate(R_SCM_getListElement(pstate, "R")));
	PROTECT(pE=duplicate(R_SCM_getListElement(pstate, "E")));
	PROTECT(pT=duplicate(R_SCM_getListElement(pstate, "T")));

	X.M = REAL(pX);
	R.M = REAL(pR);
	E.M = REAL(pE);
	T.M = REAL(pT);

	for (i=0; i<nomiss; i++) {
	  int row=MATRIX(missIdx, i, 0)-1;
	  int col=MATRIX(missIdx, i, 1)-1;
	  double e0=-eta0[col];
	  double e1=-eta1[col];
	  double thr=thresholds[col];
	  double muCond, sig2Cond;
	  double sdMu, xMax, V, proposal, r, r1, r2;
	  R_SCM_optim_data_t data={ e0, e1, 0, 0, thr };
	  double *V11;
	  int k, l;

          R_CheckUserInterrupt();

	  data.mu=muCond=MATRIX(dataMean, row, col);
	  data.sig2=sig2Cond=dataVar[col];
	  sdMu = sqrt(-1.0 / R_SCM_target_hessian(muCond, e0, e1, sig2Cond));
	  xMax = R_SCM_Brent_fmin(muCond-3*sdMu, muCond+3*sdMu,
				  R_SCM_fmin_callback,
				  &data, pow(DOUBLE_EPS, 0.25));
	  V = -1.0 / R_SCM_target_hessian(xMax, e0, e1, sig2Cond);

	  proposal = rnorm(xMax, sqrt(V));
	  r1 = R_SCM_target_dens(proposal, e0, e1, muCond, sig2Cond, thr) /
	    R_SCM_target_dens(MATRIX(X, row, col), e0, e1, muCond, sig2Cond, thr);
	  r2 = dnorm(MATRIX(X, row, col), xMax, sqrt(V), /*log=*/ 0) /
	    dnorm(proposal, xMax, sqrt(V), /*log=*/ 0);
	  r = r1 * r2;
	  if (r >= 1 || runif(0.0, 1.0) < r) {
	    noAcc++;
	    MATRIX(X, row, col) = proposal;
	    if (topresR[col]) {
	      MATRIX(R, row, col) = proposal - muCond;
	    } else if (topresE[col]) {
	      MATRIX(E, row, kj[col]-1) = proposal -  muCond;
	    } else if (topresT[col]) {
	      MATRIX(T, row, tj[col]-1) = proposal - muCond;
	    }
	  } else {
	    noRej++;
	  }
	}	

	PutRNGstate();

	PROTECT(result=NEW_LIST(6));
	PROTECT(names=NEW_CHARACTER(6));
	
	SET_VECTOR_ELT(result, 0, pX);
	SET_VECTOR_ELT(result, 1, pR);
	SET_VECTOR_ELT(result, 2, pE);
	SET_VECTOR_ELT(result, 3, pT);
	SET_VECTOR_ELT(result, 4, ScalarInteger(noAcc));
	SET_VECTOR_ELT(result, 5, ScalarInteger(noRej));

	SET_STRING_ELT(names, 0, mkChar("X"));
	SET_STRING_ELT(names, 1, mkChar("R"));
	SET_STRING_ELT(names, 2, mkChar("E"));
	SET_STRING_ELT(names, 3, mkChar("T"));
	SET_STRING_ELT(names, 4, mkChar("noAccImp"));
	SET_STRING_ELT(names, 5, mkChar("noRejImp"));

	SET_NAMES(result, names);
	
	UNPROTECT(6);
  return result;
}
