
#ifndef RSCM_H
#define RSCM_H

double R_SCM_Brent_fmin(double ax, double bx, double (*f)(double, void *),
			void *info, double tol);

#endif

