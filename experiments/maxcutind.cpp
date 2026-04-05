#include "bivariate_normal.hpp"
#include "arb_wrapper.hpp"
#include <tuple>
#include <vector>
#include <cassert>
#include <random>

double completeness = 0.8159779536207203, eps = 1e-13;
Arb comp = Arb(completeness - eps, completeness + eps);

Arb b_star = 1 - comp;
Arb s_critical = Arb::acos(-(1 - b_star) / (1 + b_star)) / Arb::pi(); // critical soundness

// slope and intercept of the affine function which our soundness wants to beat
// i.e., we want the soundness to be at least:
// slope * (bi + bj) / 2 + intercept
Arb slope =  -2 / (1 + b_star) / (1 + b_star) / Arb::pi() / Arb::sin(Arb::pi() * s_critical); 
Arb intercept = slope * (-b_star) + s_critical; 

size_t counter = 0;
Arb vol = 0;

Arb rd = 0.001;  // radius of a tiny neighborhood around b_star 
Arb x = 0.42;    // in which the slope of the threshold function is x

Arb threshs[] = {0, intercept - 0.5 + 0.001, intercept / 2 + 0.001, 0.5 - rd * x, 0.5, 0.5 + rd * x, intercept + 0.0095};
Arb brk_pts[] = {-1, -b_star, 0, b_star - rd, b_star, b_star + rd,1};

const size_t length = 7;

Arb get_t(Arb b, Arb thresh_slope, Arb thresh_intercept) { 
    // computes the threshold given the bias, 
    // and the slope/intercept corresponding to the current piece in the piecewise linear function
    return b * thresh_slope + thresh_intercept;
}

Arb get_rho(Arb bi, Arb bj) {
    // computes the relative bias of the configuration (bi, bj, -1 + bi + bj)
    if (bi.contains(-1) || bj.contains(-1)) {
        return Arb::join(-1, 0);
    }
    Arb tmp = (1 - bi) / (1 + bi) * (1 - bj) / (1 + bj);
    assert(!tmp.is_nan());
    Arb res = -Arb::safe_sqrt(tmp);
    return res;
}

// compute the first and second derivatives of rho with respect to the biases

Arb d_rho_d_b(Arb bi, Arb bj) { 
    Arb res = - 1 / Arb::safe_sqrt( -Arb::pow(1 + bi, 4) + 2 * Arb::pow(1 + bi, 3) );
    res = -res * Arb::safe_sqrt((1 - bj) / (1 + bj));
    return res;
}

Arb dd_rho_d_b2(Arb bi, Arb bj) { 
    Arb res = (-2 * Arb::pow(1 + bi, 3) + 3 * Arb::pow(1 + bi, 2)) / Arb::pow( -Arb::pow(1 + bi, 4) + 2 * Arb::pow(1 + bi, 3), Arb(3)/Arb(2));
    res = -res * Arb::safe_sqrt((1 - bj) / (1 + bj));
    return res;
}

Arb dd_rho_d_bi_d_bj(Arb bi, Arb bj) { 
    Arb res = -1 / Arb::safe_sqrt( -Arb::pow(1 + bi, 4) + 2 * Arb::pow(1 + bi, 3) );
    res = res * (-1 / Arb::safe_sqrt( -Arb::pow(1 + bj, 4) + 2 * Arb::pow(1 + bj, 3) ));
    res = -res;
    return res;
}

// -------------------------

Arb phi_rho(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) {
    // computes \phi_\rho( \Phi^{-1}(t(bi)), \Phi^{-1}(t(bj)) )

    Arb rho = get_rho(bi, bj);
    Arb ti_inv = get_t(bi, thresh_slope_i, thresh_intercept_i).norm_cdf_inv();
    Arb tj_inv = get_t(bj, thresh_slope_j, thresh_intercept_j).norm_cdf_inv();
    Arb exponent = - (ti_inv.sqr() - 2 * rho * ti_inv * tj_inv + tj_inv.sqr()) / 2 / (1 - rho.sqr());
    return 1 / Arb::pi() / 2 / Arb::safe_sqrt(1 - rho.sqr()) * Arb::exp(exponent);
}


Arb sound_actual(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) {
    // computes the soundness achieved by our rounding scheme on (bi, bj)
    Arb ti = get_t(bi, thresh_slope_i, thresh_intercept_i);
    Arb tj = get_t(bj, thresh_slope_j, thresh_intercept_j);
    Arb rho = get_rho(bi, bj);
    Arb res = 1 - biv_norm_cdf_norm_thresh(ti, tj, rho) - biv_norm_cdf_norm_thresh(1 - ti, 1 - tj, rho);
    return res;
}

// --------------------------------------------

Arb d_t_inv_d_b(Arb b, Arb thresh_slope, Arb thresh_intercept) {
    Arb t = get_t(b, thresh_slope, thresh_intercept);
    return thresh_slope / Arb::norm_pdf(t.norm_cdf_inv());
}

Arb d_phi_rho_d_bi(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) {
    Arb ti = get_t(bi, thresh_slope_i, thresh_intercept_i);
    Arb tj = get_t(bj, thresh_slope_j, thresh_intercept_j);
    Arb ti_inv = ti.norm_cdf_inv();
    Arb tj_inv = tj.norm_cdf_inv();
    Arb rho = get_rho(bi, bj);
    Arb phi_rho_res = phi_rho(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);

    Arb res = 0;

    res = d_rho_d_b(bi, bj) * 
        (rho / (1 - rho.sqr()) * phi_rho_res + 
            (ti_inv * tj_inv * (1 - rho.sqr()) - 
            rho * (ti_inv.sqr() - 2 * rho * ti_inv * tj_inv + tj_inv.sqr())) / (1 - rho.sqr()).sqr() * phi_rho_res);
    res = res + d_t_inv_d_b(bi, thresh_slope_i, thresh_intercept_i) * phi_rho_res * (-(ti_inv - rho * tj_inv) / (1 - rho.sqr()));
    return res;
}

Arb d_s_d_bi(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) { 

    Arb ti = get_t(bi, thresh_slope_i, thresh_intercept_i);
    Arb tj = get_t(bj, thresh_slope_j, thresh_intercept_j);
    Arb ti_inv = ti.norm_cdf_inv();
    Arb tj_inv = tj.norm_cdf_inv();
    Arb rho = get_rho(bi, bj);
    Arb phi_rho_res = phi_rho(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
    Arb res = 0;

    res = thresh_slope_i * ((tj_inv - rho * ti_inv) / Arb::safe_sqrt(1 - rho * rho)).norm_cdf() + phi_rho_res * d_rho_d_b(bi, bj); 
    res = (-2) * res;
    res = res + thresh_slope_i;
    
    return res;
}

Arb dd_s_d_bi2(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) {

    Arb ti = get_t(bi, thresh_slope_i, thresh_intercept_i);
    Arb tj = get_t(bj, thresh_slope_j, thresh_intercept_j);
    Arb ti_inv = ti.norm_cdf_inv();
    Arb tj_inv = tj.norm_cdf_inv();
    Arb rho = get_rho(bi, bj);
    Arb phi_rho_res = phi_rho(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
    Arb res = 0;
    Arb tmp = (tj_inv - rho * ti_inv) / Arb::safe_sqrt(1 - rho.sqr());

    res = dd_rho_d_b2(bi, bj) * phi_rho_res + d_rho_d_b(bi, bj) * d_phi_rho_d_bi(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
    res = res + thresh_slope_i * (tmp.norm_pdf()) * ((-rho * d_t_inv_d_b(bi, thresh_slope_i, thresh_intercept_i) - d_rho_d_b(bi, bj) * ti_inv) / Arb::safe_sqrt(1 - rho.sqr()) + tmp / (1 - rho.sqr()) * rho * d_rho_d_b(bi, bj));

    res = (-2) * res;

    return res;
}

Arb dd_s_d_bi_d_bj(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) {

    Arb ti = get_t(bi, thresh_slope_i, thresh_intercept_i);
    Arb tj = get_t(bj, thresh_slope_j, thresh_intercept_j);
    Arb ti_inv = ti.norm_cdf_inv();
    Arb tj_inv = tj.norm_cdf_inv();
    Arb rho = get_rho(bi, bj);
    Arb tmp = (tj_inv - rho * ti_inv) / Arb::safe_sqrt(1 - rho.sqr());
    Arb phi_rho_res = phi_rho(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);


    Arb res = 0;

    res = dd_rho_d_bi_d_bj(bi, bj) * phi_rho_res + d_rho_d_b(bi, bj) * d_phi_rho_d_bi(bj, bi, thresh_slope_j, thresh_slope_i, thresh_intercept_j, thresh_intercept_i);
    res = res + thresh_slope_i * (tmp.norm_pdf()) * ((d_t_inv_d_b(bj, thresh_slope_j, thresh_intercept_j)  - d_rho_d_b(bj, bi) * ti_inv) / Arb::safe_sqrt(1 - rho.sqr()) + tmp / (1 - rho.sqr()) * rho * d_rho_d_b(bj, bi));

    res = (-2) * res;

    return res;
}


bool recurse(Arb bi, Arb bj, Arb thresh_slope_i, Arb thresh_slope_j, Arb thresh_intercept_i, Arb thresh_intercept_j) {
    counter += 1;

    // printing volume of already checked regions
    // as well as the intervals for current biases
    printf("HERE %zu %.10f\n", counter, vol.mid_to_double());
    printf("bi = (%.10f, %.10f)\n", bi.left_edge_to_double(), bi.right_edge_to_double());
    printf("bj = (%.10f, %.10f)\n", bj.left_edge_to_double(), bj.right_edge_to_double());


    if ((bi + bj < 0) || (bi > bj)) { 
        // either infeasible or can be eliminated using symmetry
        vol = vol + bi.rad() * bj.rad();
        return true;
    }
    if ((Arb::abs(bi - b_star) <= rd) && (Arb::abs(bj - b_star) <= rd)) {
        // close enough to critical point, check Hessian
        Arb dbi2 = dd_s_d_bi2(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
        Arb dbj2 = dd_s_d_bi2(bj, bi, thresh_slope_j, thresh_slope_i, thresh_intercept_j, thresh_intercept_i);
        Arb dbidbj = dd_s_d_bi_d_bj(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
        Arb Hessian_det = dbi2 * dbj2 - dbidbj.sqr();
        if ((Hessian_det > 0) && (dbi2 > 0)) {
            printf("Hessian Good!\n");
            vol = vol + bi.rad() * bj.rad();
            return true;
        }
        else if ((Hessian_det <= 0) || (dbi2 <= 0)) {
            printf("Hessian Bad!\n");
            return false;
        }
    } 

    // compute soundness and compare it to the target
    Arb s_desired = slope * (bi + bj) / 2 + intercept;
    Arb s_actual = sound_actual(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
    printf("Soundness wanted = (%.10f, %.10f)\n", s_desired.left_edge_to_double(), s_desired.right_edge_to_double());
    printf("Soundness actual = (%.10f, %.10f)\n", s_actual.left_edge_to_double(), s_actual.right_edge_to_double());
    if (s_desired <= s_actual) {
        vol = vol + bi.rad() * bj.rad();
        return true;
    }
    else if (s_desired > s_actual) {
        printf("Counterexample found!\n");
        vol = vol + bi.rad() * bj.rad();
        return false;
    }
    else {
        // since we are comparing intervals, it could be that neither of the previous two cases holds

        if (s_actual.rad() < 1e-10) {
            // error bar for soundness is too small
            // stop recursion and return inconclusive
            printf("Suspicious point found!\n");
            return false;
        }
        if (bi.rad() > 1e-7) {
            // check derivative in b_i and optimize if it is positive/negative
            Arb dbi = d_s_d_bi(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j) - slope / 2;
            if (dbi > 0) {
                vol = vol + bi.rad() * bj.rad();
                bi = bi.left_edge();
            } else if (dbi < 0) {
                vol = vol + bi.rad() * bj.rad();
                bi = bi.right_edge();
            }
        }

        if (bj.rad() > 1e-7) {      
            // check derivative in b_j and optimize if it is positive/negative 
            Arb dbj = d_s_d_bi(bj, bi, thresh_slope_j, thresh_slope_i, thresh_intercept_j, thresh_intercept_i) - slope / 2;
            if (dbj > 0) {
                vol = vol + bi.rad() * bj.rad();
                bj = bj.left_edge();
            } else if (dbj < 0) {
                vol = vol + bi.rad() * bj.rad();
                bj = bj.right_edge();
            }
        }

        // recurse
        Arb bi_l = bi.left_half();
        Arb bi_r = bi.right_half();
        Arb bj_l = bj.left_half();
        Arb bj_r = bj.right_half();
        return recurse(bi_l, bj_l, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j) && 
               recurse(bi_l, bj_r, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j) && 
               recurse(bi_r, bj_l, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j) && 
               recurse(bi_r, bj_r, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
    }
    
    return true;
}

int main(int argc, char* argv[]) {
    comp.println();
    
    Arb b_limit = Arb(-1,1);
    Arb bi;
    Arb bj;

    Arb thresh_slope_i;
    Arb thresh_slope_j;

    Arb thresh_intercept_i;
    Arb thresh_intercept_j;

    bool res = true;

    for (size_t i = 0; i < length - 1; i++) {
        bi = Arb::join(brk_pts[i], brk_pts[i+1]);
        bi = bi.intersect(b_limit);

        thresh_slope_i = (threshs[i+1] - threshs[i]) / (brk_pts[i + 1] - brk_pts[i]);
        thresh_intercept_i = threshs[i] - thresh_slope_i * brk_pts[i];


        for (size_t j = i; j < length - 1; j++) {
            bj = Arb::join(brk_pts[j], brk_pts[j+1]);
            bj = bj.intersect(b_limit);

            thresh_slope_j = (threshs[j+1] - threshs[j]) / (brk_pts[j + 1] - brk_pts[j]);
            thresh_intercept_j = threshs[j] - thresh_slope_j * brk_pts[j];
            res = recurse(bi, bj, thresh_slope_i, thresh_slope_j, thresh_intercept_i, thresh_intercept_j);
            if (!res)
                break;
        }
        if (!res)
            break;
    }

    printf("ANSWER %d\n", res);

    return 0;
}