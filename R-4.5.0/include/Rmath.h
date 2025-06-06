/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2025  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.
   Some entries possibly are not yet documented in 'Writing R Extensions'.

*/
#ifndef RMATH_H
#define RMATH_H

/* needed for cospi etc */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
# include <cmath>
#else
# include <math.h>
#endif

#ifdef NO_C_HEADERS
# warning "use of NO_C_HEADERS is defunct and will be ignored"
#endif

/*-- Mathlib as part of R --  define this for standalone : */
/* #undef MATHLIB_STANDALONE */

#define R_VERSION_STRING "4.5.0"

// Legacy defines -- C99 functions which R >= 3.5.0 requires
#ifndef HAVE_EXPM1
# define HAVE_EXPM1 1
#endif
#ifndef HAVE_HYPOT
# define HAVE_HYPOT 1
#endif
#ifndef HAVE_LOG1P
# define HAVE_LOG1P 1
#endif

#ifndef HAVE_WORKING_LOG1P
# define HAVE_WORKING_LOG1P 1
#endif

#if !defined(HAVE_WORKING_LOG1P)
/* remap to avoid problems with getting the right entry point */
extern "C" {
double  Rlog1p(double);
}
#define log1p Rlog1p
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


#ifdef MATHLIB_STANDALONE
# ifndef R_EXT_BOOLEAN_H_
/* "copy-paste" R_ext/Boolean.h if not already included: 
   This is standalone, so we do not worry about the size
*/
 #define R_EXT_BOOLEAN_H_
 #undef FALSE
 #undef TRUE
 typedef enum { FALSE = 0, TRUE } Rboolean;
 // stdbool.h is obsolescent but we could pro tem just use it.
 // But to future-proof copy from Boolean.h
 #if defined __STDC_VERSION__ && __STDC_VERSION__ > 202000L
 #elif defined __cplusplus
 #else
  # include <stdbool.h>
 #endif

# endif
#else
# include <R_ext/Boolean.h>
#endif

#if !defined(MATHLIB_STANDALONE)
#define bessel_i	Rf_bessel_i
#define bessel_j	Rf_bessel_j
#define bessel_k	Rf_bessel_k
#define bessel_y	Rf_bessel_y
#define bessel_i_ex	Rf_bessel_i_ex
#define bessel_j_ex	Rf_bessel_j_ex
#define bessel_k_ex	Rf_bessel_k_ex
#define bessel_y_ex	Rf_bessel_y_ex
#define beta		Rf_beta
#define choose		Rf_choose
#define dbeta		Rf_dbeta
#define dbinom		Rf_dbinom
#define dbinom_raw	Rf_dbinom_raw
#define dcauchy		Rf_dcauchy
#define dchisq		Rf_dchisq
#define dexp		Rf_dexp
#define df		Rf_df
#define dgamma		Rf_dgamma
#define dgeom		Rf_dgeom
#define dhyper		Rf_dhyper
#define digamma		Rf_digamma
#define dlnorm		Rf_dlnorm
#define dlogis		Rf_dlogis
#define dnbeta		Rf_dnbeta
#define dnbinom		Rf_dnbinom
#define dnbinom_mu	Rf_dnbinom_mu
#define dnchisq		Rf_dnchisq
#define dnf		Rf_dnf
#define dnorm4		Rf_dnorm4
#define dnt		Rf_dnt
#define dpois_raw	Rf_dpois_raw
#define dpois		Rf_dpois
#define dpsifn		Rf_dpsifn
#define dsignrank	Rf_dsignrank
#define dt		Rf_dt
#define dtukey		Rf_dtukey
#define dunif		Rf_dunif
#define dweibull	Rf_dweibull
#define dwilcox		Rf_dwilcox
#define fmax2		Rf_fmax2
#define fmin2		Rf_fmin2
#define fprec		Rf_fprec
#define fround		Rf_fround
#define ftrunc		Rf_ftrunc
#define fsign		Rf_fsign
#define gammafn		Rf_gammafn
#define imax2		Rf_imax2
#define imin2		Rf_imin2
#define lbeta		Rf_lbeta
#define lchoose		Rf_lchoose
#define lgammafn	Rf_lgammafn
#define lgammafn_sign	Rf_lgammafn_sign
#define lgamma1p	Rf_lgamma1p
#define pow1p		Rf_pow1p
#define log1mexp       	Rf_log1mexp
#define log1pexp       	Rf_log1pexp
#define log1pmx		Rf_log1pmx
#define logspace_add	Rf_logspace_add
#define logspace_sub	Rf_logspace_sub
#define logspace_sum	Rf_logspace_sum
#define pbeta		Rf_pbeta
#define pbeta_raw	Rf_pbeta_raw
#define pbinom		Rf_pbinom
#define pcauchy		Rf_pcauchy
#define pchisq		Rf_pchisq
#define pentagamma	Rf_pentagamma
#define pexp		Rf_pexp
#define pf		Rf_pf
#define pgamma		Rf_pgamma
#define pgeom		Rf_pgeom
#define phyper		Rf_phyper
#define plnorm		Rf_plnorm
#define plogis		Rf_plogis
#define pnbeta		Rf_pnbeta
#define pnbinom		Rf_pnbinom
#define pnbinom_mu     	Rf_pnbinom_mu
#define pnchisq		Rf_pnchisq
#define pnf		Rf_pnf
#define pnorm5		Rf_pnorm5
#define pnorm_both	Rf_pnorm_both
#define pnt		Rf_pnt
#define ppois		Rf_ppois
#define psignrank	Rf_psignrank
#define psigamma	Rf_psigamma
#define pt		Rf_pt
#define ptukey		Rf_ptukey
#define punif		Rf_punif
#define pweibull	Rf_pweibull
#define pwilcox		Rf_pwilcox
#define qbeta		Rf_qbeta
#define qbinom		Rf_qbinom
#define qcauchy		Rf_qcauchy
#define qchisq		Rf_qchisq
#define qchisq_appr	Rf_qchisq_appr
#define qexp		Rf_qexp
#define qf		Rf_qf
#define qgamma		Rf_qgamma
#define qgeom		Rf_qgeom
#define qhyper		Rf_qhyper
#define qlnorm		Rf_qlnorm
#define qlogis		Rf_qlogis
#define qnbeta		Rf_qnbeta
#define qnbinom		Rf_qnbinom
#define qnbinom_mu     	Rf_qnbinom_mu
#define qnchisq		Rf_qnchisq
#define qnf		Rf_qnf
#define qnorm5		Rf_qnorm5
#define qnt		Rf_qnt
#define qpois		Rf_qpois
#define qsignrank	Rf_qsignrank
#define qt		Rf_qt
#define qtukey		Rf_qtukey
#define qunif		Rf_qunif
#define qweibull	Rf_qweibull
#define qwilcox		Rf_qwilcox
#define rbeta		Rf_rbeta
#define rbinom		Rf_rbinom
#define rcauchy		Rf_rcauchy
#define rchisq		Rf_rchisq
#define rexp		Rf_rexp
#define rf		Rf_rf
#define rgamma		Rf_rgamma
#define rgeom		Rf_rgeom
#define rhyper		Rf_rhyper
#define rlnorm		Rf_rlnorm
#define rlogis		Rf_rlogis
#define rmultinom	Rf_rmultinom
#define rnbeta		Rf_rnbeta
#define rnbinom		Rf_rnbinom
#define rnbinom_mu     	Rf_rnbinom_mu
#define rnchisq		Rf_rnchisq
#define rnf		Rf_rnf
#define rnorm		Rf_rnorm
#define rnt		Rf_rnt
#define rpois		Rf_rpois
#define rsignrank	Rf_rsignrank
#define rt		Rf_rt
#define rtukey		Rf_rtukey
#define runif		Rf_runif
#define rweibull	Rf_rweibull
#define rwilcox		Rf_rwilcox
#define sign		Rf_sign
#define tetragamma	Rf_tetragamma
#define trigamma	Rf_trigamma
#endif /* if not STANDALONE */

#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
extern "C" {
#endif
	/* R's versions with !R_FINITE checks */

double R_pow(double x, double y);
double R_pow_di(double, int);

	/* Random Number Generators */

double	norm_rand(void); //not remapped
double	unif_rand(void); //not remapped
double  R_unif_index(double);
double	exp_rand(void); //not remapped
#ifdef MATHLIB_STANDALONE
void	set_seed(unsigned int, unsigned int);
void	get_seed(unsigned int *, unsigned int *);
#endif

	/* Normal Distribution */

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
double	rnorm(double, double);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	dunif(double, double, double, int);
double	punif(double, double, double, int, int);
double	qunif(double, double, double, int, int);
double	runif(double, double);

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);
double	rgamma(double, double);

double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
double  log1pexp(double); // <-- ../nmath/plogis.c
double  log1mexp(double);
double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

double  pow1p(double, double); /* pow1p(x, y) := (1+x)^y  accurately also for |x| << 1 */

/* Compute the log of a sum or difference from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 * or  log (exp (logx) - exp (logy))
 *
 * without causing overflows or throwing away too much accuracy:
 */
double  logspace_add(double logx, double logy);
double  logspace_sub(double logx, double logy);
double  logspace_sum(const double *, int);

	/* Beta Distribution */

double	dbeta(double, double, double, int);
double	pbeta(double, double, double, int, int);
double	qbeta(double, double, double, int, int);
double	rbeta(double, double);

	/* Lognormal Distribution */

double	dlnorm(double, double, double, int);
double	plnorm(double, double, double, int, int);
double	qlnorm(double, double, double, int, int);
double	rlnorm(double, double);

	/* Chi-squared Distribution */

double	dchisq(double, double, int);
double	pchisq(double, double, int, int);
double	qchisq(double, double, int, int);
double	rchisq(double);

	/* Non-central Chi-squared Distribution */

double	dnchisq(double, double, double, int);
double	pnchisq(double, double, double, int, int);
double	qnchisq(double, double, double, int, int);
double	rnchisq(double, double);

	/* F Distribution */

double	df(double, double, double, int);
double	pf(double, double, double, int, int);
double	qf(double, double, double, int, int);
double	rf(double, double);

	/* Student t Distribution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);
double	rt(double);

	/* Binomial Distribution */

double  dbinom_raw(double x, double n, double p, double q, int give_log);
double	dbinom(double, double, double, int);
double	pbinom(double, double, double, int, int);
double	qbinom(double, double, double, int, int);
double	rbinom(double, double);

	/* Multinomial Distribution */

void	rmultinom(int, double*, int, int*);

	/* Cauchy Distribution */

double	dcauchy(double, double, double, int);
double	pcauchy(double, double, double, int, int);
double	qcauchy(double, double, double, int, int);
double	rcauchy(double, double);

	/* Exponential Distribution */

double	dexp(double, double, int);
double	pexp(double, double, int, int);
double	qexp(double, double, int, int);
double	rexp(double);

	/* Geometric Distribution */

double	dgeom(double, double, int);
double	pgeom(double, double, int, int);
double	qgeom(double, double, int, int);
double	rgeom(double);

	/* Hypergeometric Distribution */

double	dhyper(double, double, double, double, int);
double	phyper(double, double, double, double, int, int);
double	qhyper(double, double, double, double, int, int);
double	rhyper(double, double, double);

	/* Negative Binomial Distribution */

double	dnbinom(double, double, double, int);
double	pnbinom(double, double, double, int, int);
double	qnbinom(double, double, double, int, int);
double	rnbinom(double, double);

double	dnbinom_mu(double, double, double, int);
double	pnbinom_mu(double, double, double, int, int);
double	qnbinom_mu(double, double, double, int, int);
double	rnbinom_mu(double, double);

	/* Poisson Distribution */

double	dpois_raw (double, double, int);
double	dpois(double, double, int);
double	ppois(double, double, int, int);
double	qpois(double, double, int, int);
double	rpois(double);

	/* Weibull Distribution */

double	dweibull(double, double, double, int);
double	pweibull(double, double, double, int, int);
double	qweibull(double, double, double, int, int);
double	rweibull(double, double);

	/* Logistic Distribution */

double	dlogis(double, double, double, int);
double	plogis(double, double, double, int, int);
double	qlogis(double, double, double, int, int);
double	rlogis(double, double);

	/* Non-central Beta Distribution */

double	dnbeta(double, double, double, double, int);
double	pnbeta(double, double, double, double, int, int);
double	qnbeta(double, double, double, double, int, int);
double	rnbeta(double, double, double);

	/* Non-central F Distribution */

double  dnf(double, double, double, double, int);
double	pnf(double, double, double, double, int, int);
double	qnf(double, double, double, double, int, int);

	/* Non-central Student t Distribution */

double	dnt(double, double, double, int);
double	pnt(double, double, double, int, int);
double	qnt(double, double, double, int, int);

	/* Studentized Range Distribution */

double	ptukey(double, double, double, double, int, int);
double	qtukey(double, double, double, double, int, int);

	/* Wilcoxon Rank Sum Distribution */

double dwilcox(double, double, double, int);
double pwilcox(double, double, double, int, int);
double qwilcox(double, double, double, int, int);
double rwilcox(double, double);
void wilcox_free(void); // not remapped
	/* Wilcoxon Signed Rank Distribution */

double dsignrank(double, double, int);
double psignrank(double, double, int, int);
double qsignrank(double, double, int, int);
double rsignrank(double);
void signrank_free(void); // not remapped

	/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);
void    dpsifn(double, int, int, int, double*, int*, int*);
double	psigamma(double, double);
double	digamma(double);
double	trigamma(double);
double	tetragamma(double);
double	pentagamma(double);

double	beta(double, double);
double	lbeta(double, double);

double	choose(double, double);
double	lchoose(double, double);

	/* Bessel Functions */

double	bessel_i(double, double, double);
double	bessel_j(double, double);
double	bessel_k(double, double, double);
double	bessel_y(double, double);
double	bessel_i_ex(double, double, double, double *);
double	bessel_j_ex(double, double, double *);
double	bessel_k_ex(double, double, double, double *);
double	bessel_y_ex(double, double, double *);


	/* General Support Functions */

int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
double	sign(double);
double	fprec(double, double);
double	fround(double, double);
double	fsign(double, double);
double	ftrunc(double);

/* More accurate cos(pi*x), sin(pi*x), tan(pi*x)

   These declarations might clash with system headers if someone had
   already included math.h with __STDC_WANT_IEC_60559_FUNCS_EXT__
   defined (and we try, above).
   We check for that via the value of __STDC_IEC_60559_FUNCS__
*/
#if !(defined(__STDC_IEC_60559_FUNCS__) && __STDC_IEC_60559_FUNCS__ >= 201506L)
double cospi(double);
double sinpi(double);
double tanpi(double);
#endif
double Rtanpi(double); /* our own in any case */

#if !defined(MATHLIB_STANDALONE) && defined(R_NO_REMAP_RMATH)
#undef bessel_i
#undef bessel_j
#undef bessel_k
#undef bessel_y
#undef bessel_i_ex
#undef bessel_j_ex
#undef bessel_k_ex
#undef bessel_y_ex
#undef beta
#undef choose
#undef dbeta
#undef dbinom
#undef dbinom_raw
#undef dcauchy
#undef dchisq
#undef dexp
#undef df
#undef dgamma
#undef dgeom
#undef dhyper
#undef digamma
#undef dlnorm
#undef dlogis
#undef dnbeta
#undef dnbinom
#undef dnbinom_mu
#undef dnchisq
#undef dnf
#undef dnorm4
#undef dnt
#undef dpois_raw
#undef dpois
#undef dpsifn
#undef dsignrank
#undef dt
#undef dtukey
#undef dunif
#undef dweibull
#undef dwilcox
#undef fmax2
#undef fmin2
#undef fprec
#undef fround
#undef ftrunc
#undef fsign
#undef gammafn
#undef imax2
#undef imin2
#undef lbeta
#undef lchoose
#undef lgammafn
#undef lgammafn_sign
#undef lgamma1p
#undef pow1p
#undef log1mexp
#undef log1pexp
#undef log1pmx
#undef logspace_add
#undef logspace_sub
#undef logspace_sum
#undef pbeta
#undef pbeta_raw
#undef pbinom
#undef pcauchy
#undef pchisq
#undef pentagamma
#undef pexp
#undef pf
#undef pgamma
#undef pgeom
#undef phyper
#undef plnorm
#undef plogis
#undef pnbeta
#undef pnbinom
#undef pnbinom_mu
#undef pnchisq
#undef pnf
#undef pnorm5
#undef pnorm_both
#undef pnt
#undef ppois
#undef psignrank
#undef psigamma
#undef pt
#undef ptukey
#undef punif
#undef pweibull
#undef pwilcox
#undef qbeta
#undef qbinom
#undef qcauchy
#undef qchisq
#undef qchisq_appr
#undef qexp
#undef qf
#undef qgamma
#undef qgeom
#undef qhyper
#undef qlnorm
#undef qlogis
#undef qnbeta
#undef qnbinom
#undef qnbinom_mu
#undef qnchisq
#undef qnf
#undef qnorm5
#undef qnt
#undef qpois
#undef qsignrank
#undef qt
#undef qtukey
#undef qunif
#undef qweibull
#undef qwilcox
#undef rbeta
#undef rbinom
#undef rcauchy
#undef rchisq
#undef rexp
#undef rf
#undef rgamma
#undef rgeom
#undef rhyper
#undef rlnorm
#undef rlogis
#undef rmultinom
#undef rnbeta
#undef rnbinom
#undef rnbinom_mu
#undef rnchisq
#undef rnf
#undef rnorm
#undef rnt
#undef rpois
#undef rsignrank
#undef rt
#undef rtukey
#undef runif
#undef rweibull
#undef rwilcox
#undef sign
#undef tetragamma
#undef trigamma

#undef dnorm
#undef pnorm
#undef qnorm
#endif /* not STANDALONE *and* R_NO_REMAP_RMATH */

/* ----------------- Private part of the header file ------------------- */

#if defined(MATHLIB_STANDALONE) && !defined(MATHLIB_PRIVATE_H)
/* second is defined by nmath.h */

/* If isnan is a macro, as C99 specifies, the C++
   math header will undefine it. This happens on macOS */
# ifdef __cplusplus
  int R_isnancpp(double); /* in mlutils.c */
#  define ISNAN(x)     R_isnancpp(x)
# else
#  define ISNAN(x)     (isnan(x)!=0)
# endif

# define R_FINITE(x)    R_finite(x)
int R_finite(double);

# ifdef _WIN32  /* not Win32 as no config information */
#  ifdef RMATH_DLL
#   define R_EXTERN extern __declspec(dllimport)
#  else
#   define R_EXTERN extern
#  endif
R_EXTERN double NA_REAL;
R_EXTERN double R_PosInf;
R_EXTERN double R_NegInf;
R_EXTERN int N01_kind;
#  undef R_EXTERN
#else
extern int N01_kind;
# endif

#endif /* MATHLIB_STANDALONE and not PRIVATE_H */

#ifdef  __cplusplus
}
#endif

#endif /* RMATH_H */
