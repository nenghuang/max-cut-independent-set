#include "bivariate_normal.hpp"
#include "flint/acb_calc.h"
#include <cassert>

Arb biv_norm_cdf_norm_thresh(const Arb &t1, const Arb &t2, const Arb &rho) {
    // key observation: the problem is monotone in all three parameters
    assert(!t1.is_nan());
    assert(!t2.is_nan());
    assert(!rho.is_nan());

    Arb lower, upper;

    if (t1.right_edge() >= T_UPPER_THRESH) {
        upper = t2.right_edge();
    }
    else if (t2.right_edge() >= T_UPPER_THRESH) {
        upper = t1.right_edge();
    }
    else {
        Arb t1u = Arb::max(t1.right_edge(), T_LOWER_THRESH);
        Arb t2u = Arb::max(t2.right_edge(), T_LOWER_THRESH);
        upper = biv_norm_cdf(t1u.norm_cdf_inv(), t2u.norm_cdf_inv(), rho);
    }

    if (t1.left_edge() <= T_LOWER_THRESH ||
        t2.left_edge() <= T_LOWER_THRESH) {
        lower = 0;
    }
    else {
        Arb t1l = Arb::min(t1.left_edge(), T_UPPER_THRESH);
        Arb t2l = Arb::min(t2.left_edge(), T_UPPER_THRESH);
        lower = biv_norm_cdf(t1l.norm_cdf_inv(), t2l.norm_cdf_inv(), rho);
    }

    return Arb::join(lower, upper);
}


// any part of rho extending beyond +/-1 is ignored
Arb biv_norm_cdf(const Arb &t1, const Arb &t2, const Arb &rho) {
    // key observation: the problem is monotone in all three parameters

    /* flint_printf("BIV NORM DEBUG\n");
    flint_printf("t1: "); t1.println();
    flint_printf("t2: "); t2.println();
    flint_printf("rho: "); rho.println(); */
 
    assert(!t1.is_nan());
    assert(!t2.is_nan());
    assert(!rho.is_nan());

    Arb lower, upper;

    // snap upper/lower end to +/- 1
    if (rho.right_edge() >= RHO_THRESH) {
        upper = biv_norm_cdf_rho_plus_one(t1.right_edge(), t2.right_edge());
    }
    else if (rho.right_edge() <= -RHO_THRESH) {
        upper = biv_norm_cdf_unsafe(t1.right_edge(), t2.right_edge(), -RHO_THRESH);
    }
    else {
        upper = biv_norm_cdf_unsafe(t1.right_edge(), t2.right_edge(), rho.right_edge());
    }

    if (rho.left_edge() <= -RHO_THRESH) {
        lower = biv_norm_cdf_rho_minus_one(t1.left_edge(), t2.left_edge());
    }
    else if (rho.left_edge() >= RHO_THRESH) {
        lower = biv_norm_cdf_unsafe(t1.left_edge(), t2.left_edge(), RHO_THRESH);
    }
    else {
        lower = biv_norm_cdf_unsafe(t1.left_edge(), t2.left_edge(), rho.left_edge());
    }
    
    return Arb::join(lower, upper);
}

Arb biv_norm_cdf_rho_plus_one(const Arb &t1, const Arb &t2) {
    // min(cdf(t1), cdf(t2))
    return Arb::min(t1.norm_cdf(), t2.norm_cdf());
}

Arb biv_norm_cdf_rho_minus_one(const Arb &t1, const Arb &t2) {
    // max(0, 1 - cdf(-t1) - cdf(-t2)) 
    // = max(0, cdf(t1) + cdf(t2) -1 )
    return Arb::max(0, t1.norm_cdf() + t2.norm_cdf() - 1);
}

Arb biv_norm_cdf_unsafe(const Arb &t1, const Arb &t2, const Arb &rho) {
    // special handling when rho is close to +/-1 is NOT implemented...

    acb_calc_integrate_opt_t options;
    acb_calc_integrate_opt_init(options);

    Acb ct1(t1), ct2(t2);
    acb_t param[2] = {(*ct1.t), (*ct2.t)};

    mag_t tol; mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -GLOBAL_PRECISION);

    slong goal = GLOBAL_PRECISION;
    slong prec = GLOBAL_PRECISION;

    Acb a(0), b(rho), res;

    acb_calc_integrate(res.t, _biv_norm_cdf_helper, param, a.t, b.t, goal, tol, options, prec);

    assert(res.is_real());

    Arb ans = res.real();
    
    // Drezner & Wesolowsky eq.6 with additive correction
    return ans + t1.norm_cdf() * t2.norm_cdf(); 
}

Arb biv_norm_cdf_d_t1(const Arb &t1, const Arb &t2, const Arb &rho) {
    // formula is norm_pdf(t1) * norm_cdf((t2 - rho * t1) / sqrt (1 - rho*rho))
    Arb u = t1.norm_pdf();
    Arb v = (t2 - rho * t1) / Arb::sqrt(1 - rho.sqr());
    return u*v.norm_cdf(); 
}

Arb biv_norm_cdf_d_t2(const Arb &t1, const Arb &t2, const Arb &rho) {
    return biv_norm_cdf_d_t1(t2, t1, rho); 
}

Arb biv_norm_cdf_d_rho(const Arb &t1, const Arb &t2, const Arb &rho) {
    // 1/(2*pi*sqrt(1-rho*rho)) * exp(- (t1*t1 - 2*rho*t1*t2 + t2*t2)/ (2*(1-rho*rho)))
    Arb a = 1 - rho.sqr();
    Arb b = t1.sqr() + t2.sqr() - 2*rho*t1*t2;
    return Arb::exp((-1.0/2.0)*b/a) / (2*Arb::pi()*a.sqrt());
}

int _biv_norm_cdf_helper(acb_ptr res, const acb_t rho, void * param, slong order, slong prec) {
    // documentation says these should never be tripped...
    assert(order == 0 || order == 1);
    assert(prec <= GLOBAL_PRECISION);

    acb_t* param_acb_t = (acb_t*) param;
    Acb t1(param_acb_t[0]), t2(param_acb_t[1]);

    Acb ans = biv_norm_cdf_d_rho_analytic(t1, t2, rho, order != 0);

    acb_set(res, ans.t);

    return 0;
}

Acb biv_norm_cdf_d_rho_analytic(const Acb &t1, const Acb &t2, const Acb &rho, int analytic) {
    // 1/(2*pi*sqrt(1-rho*rho)) * exp(- (t1*t1 - 2*rho*t1*t2 + t2*t2)/ (2*(1-rho*rho)))
    Acb a = 1 - rho.sqr();
    Acb b = t1.sqr() + t2.sqr() - 2*rho*t1*t2;
    return Acb::exp((-1.0/2.0)*b/a) / (2*Acb::pi()*a.sqrt_analytic(analytic));
}
