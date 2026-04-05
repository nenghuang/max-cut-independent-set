#ifndef BIVARIATE_NORMAL_HPP
#define BIVARIATE_NORMAL_HPP

#include "arb_wrapper.hpp"
#include "acb_wrapper.hpp"

// given a bivariate normal with covariance matrix (1 rho; rho 1)
// find the probability that the sample (X,Y) satisfies
// X <= t1 and Y <= t2
// ``unsafe'' version has improper handling when rho is +/- 1
// this parameter decides when we should switch to special treatment
const Arb RHO_THRESH(1-1e-16);
const Arb T_UPPER_THRESH(1-1e-16);
const Arb T_LOWER_THRESH(1e-16);

Arb biv_norm_cdf_norm_thresh(const Arb &t1, const Arb &t2, const Arb &rho);

Arb biv_norm_cdf(const Arb &t1, const Arb &t2, const Arb &rho);
Arb biv_norm_cdf_unsafe(const Arb &t1, const Arb &t2, const Arb &rho);
Arb biv_norm_cdf_rho_plus_one(const Arb &t1, const Arb &t2);
Arb biv_norm_cdf_rho_minus_one(const Arb &t1, const Arb &t2);

// partial derivatives
Arb biv_norm_cdf_d_t1(const Arb &t1, const Arb &t2, const Arb &rho);
Arb biv_norm_cdf_d_t2(const Arb &t1, const Arb &t2, const Arb &rho);
Arb biv_norm_cdf_d_rho(const Arb &t1, const Arb &t2, const Arb &rho);

// might add second order partials later

// helper function for biv_norm_cdf
int _biv_norm_cdf_helper(acb_ptr res, const acb_t rho, void * param, slong order, slong prec);

// need Acb analytic version, too!
Acb biv_norm_cdf_d_rho_analytic(const Acb &t1, const Acb &t2, const Acb &rho, int analytic);
#endif
