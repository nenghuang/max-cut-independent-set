#include "arb_wrapper.hpp"
#include "flint/arf.h"
#include "flint/arb_hypgeom.h"
#include <cassert>

// constructors

Arb::Arb() {
    arb_init(this -> t);
}

Arb::Arb(double d) {
    arb_init(this -> t);
    arb_set_d(this -> t, d);
}

Arb::Arb(double d1, double d2) {
    arb_init(this -> t);
    Arb x(d1), y(d2);
    arb_union(this -> t, x.t, y.t, GLOBAL_PRECISION);
}

// destructor

Arb::~Arb() {
    arb_clear(this -> t);
}

// display info

void Arb::print() const {
    arb_printd(this -> t, 10);
}

void Arb::println() const {
    this -> print();
    flint_printf("\n");
}

void Arb::pretty_print() const {
    flint_printf("[");
    arb_printd((this -> left_edge()).t, 10);
    flint_printf(", ");
    arb_printd((this -> right_edge()).t, 10);
    flint_printf("]");
}

void Arb::pretty_println() const {
    this -> pretty_print();
    flint_printf("\n");
}

// export
double Arb::mid_to_double() const {
    arf_t u;
    arf_init(u);
    arb_get_ubound_arf(u, (this->mid()).t, GLOBAL_PRECISION);
    double ans = arf_get_d(u, ARF_RND_NEAR);
    arf_clear(u);
    return ans;
}

double Arb::left_edge_to_double() const {
    arf_t u;
    arf_init(u);
    arb_get_ubound_arf(u, (this->left_edge()).t, GLOBAL_PRECISION);
    double ans = arf_get_d(u, ARF_RND_NEAR);
    arf_clear(u);
    return ans;
}

double Arb::right_edge_to_double() const {
    arf_t u;
    arf_init(u);
    arb_get_ubound_arf(u, (this->right_edge()).t, GLOBAL_PRECISION);
    double ans = arf_get_d(u, ARF_RND_NEAR);
    arf_clear(u);
    return ans;
}


// interval arithmetic operations
Arb Arb::join(const Arb& rhs) const {
    Arb ans;
    arb_union(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::join(const Arb& lhs, const Arb& rhs) {
    return lhs.join(rhs);
}

Arb Arb::intersect(const Arb& rhs) const {
    Arb ans;
    if (arb_intersection(ans.t, this->t, rhs.t, GLOBAL_PRECISION)) {
        return ans;
    }
    else {
        return nan();
    }
}

Arb Arb::intersect(const Arb& lhs, const Arb& rhs) {
    return lhs.intersect(rhs);
}

int Arb::contains(const Arb& rhs) const {
    return arb_contains(this -> t, rhs.t);
}

int Arb::contains(const Arb &lhs, const Arb& rhs){
    return lhs.contains(rhs);
}

int Arb::inside(const Arb& rhs) const {
    return rhs.contains(*this);
}

int Arb::inside(const Arb &lhs, const Arb& rhs) {
    return rhs.contains(lhs);
}

Arb Arb::mid() const {
    Arb ans;
    arb_get_mid_arb(ans.t, this->t);
    return ans;
}

Arb Arb::rad() const {
    Arb ans;
    arb_get_rad_arb(ans.t, this->t);
    return ans;
}

Arb Arb::left_edge() const {
    Arb m = this->mid();
    Arb r = this->rad();
    return m-r;
}

Arb Arb::right_edge() const {
    Arb m = this->mid();
    Arb r = this->rad();
    return m+r;
}

Arb Arb::left_half() const {
    return join(this->left_edge(), this->mid());
}

Arb Arb::right_half() const {
    return join(this->mid(), this->right_edge());
}

Arb Arb::nan() {
    Arb ans;
    arb_indeterminate(ans.t);
    return ans;
}

int Arb::is_nan() const {
    return arf_is_nan(arb_midref(this->t));
}

// operators

Arb Arb::operator-() const {
    Arb ans;
    arb_neg(ans.t, this->t);
    return ans;
}


Arb Arb::operator+(const Arb& rhs) const {
    Arb ans;
    arb_add(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::operator-(const Arb& rhs) const {
    Arb ans;
    arb_sub(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::operator*(const Arb& rhs) const {
    Arb ans;
    arb_mul(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::operator/(const Arb& rhs) const {
    Arb ans;
    arb_div(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::operator+(double rhs) const {
    Arb x(rhs);
    return (*this) + x;
}

Arb Arb::operator-(double rhs) const {
    Arb x(rhs);
    return (*this) - x;
}

Arb Arb::operator*(double rhs) const {
    Arb x(rhs);
    return (*this) * x;
}

Arb Arb::operator/(double rhs) const {
    Arb x(rhs);
    return (*this) / x;
}


Arb operator+(double lhs, const Arb& rhs) {
    Arb x(lhs);
    return x + rhs;
}

Arb operator-(double lhs, const Arb& rhs) {
    Arb x(lhs);
    return x - rhs;
}

Arb operator*(double lhs, const Arb& rhs) {
    Arb x(lhs);
    return x * rhs;
}

Arb operator/(double lhs, const Arb& rhs) {
    Arb x(lhs);
    return x / rhs;
}



int Arb::operator==(const Arb& rhs) const {
    return arb_eq(this -> t, rhs.t);
}

int Arb::operator!=(const Arb& rhs) const {
    return arb_ne(this -> t, rhs.t);
}

int Arb::operator<(const Arb& rhs) const {
    return arb_lt(this -> t, rhs.t);
}

int Arb::operator<=(const Arb& rhs) const {
    return arb_le(this -> t, rhs.t);
}

int Arb::operator>(const Arb& rhs) const {
    return arb_gt(this -> t, rhs.t);
}

int Arb::operator>=(const Arb& rhs) const {
    return arb_ge(this -> t, rhs.t);
}

// mathematical constants and functions

Arb Arb::pi() {
    Arb ans;
    arb_const_pi(ans.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::abs() const {
    Arb ans;
    arb_abs(ans.t, this->t);
    return ans;
}

Arb Arb::abs(const Arb& x) {
    return x.abs();
}

Arb Arb::min(const Arb& rhs) const {
    Arb ans;
    arb_min(ans.t, this -> t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::min(const Arb& lhs, const Arb& rhs) {
    return lhs.min(rhs);
}

Arb Arb::max(const Arb& rhs) const {
    Arb ans;
    arb_max(ans.t, this -> t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::max(const Arb& lhs, const Arb& rhs) {
    return lhs.max(rhs);
}


Arb Arb::exp() const {
    Arb ans;
    arb_exp(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::exp(const Arb& x) {
    return x.exp();
}

Arb Arb::sqrt() const {
    Arb ans;
    arb_sqrt(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::sqrt(const Arb& x) {
    return x.sqrt();
}

Arb Arb::safe_sqrt() const {
    if (*this >= 0) {
        return this -> sqrt();
    }
    else if (this -> right_edge() >= 0) {
        return join(0, (this -> right_edge()).sqrt());
    }
    else {
        return 0;
    }
}

Arb Arb::safe_sqrt(const Arb& x) {
    return x.safe_sqrt();
}

Arb Arb::sqr() const {
    Arb ans;
    arb_sqr(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::sqr(const Arb& x) {
    return x.sqr();
}

Arb Arb::pow(const Arb& rhs) const {
    Arb ans;
    arb_pow(ans.t, this -> t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::pow(const Arb& lhs, const Arb& rhs) {
    return lhs.pow(rhs);
}

Arb Arb::erf() const {
    Arb ans;
    arb_hypgeom_erf(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::erf(const Arb& x) {
    return x.erf();
}

Arb Arb::erf_inv() const {
    Arb ans;
    arb_hypgeom_erfinv(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Arb Arb::erf_inv(const Arb& x) {
    return x.erf_inv();
}

Arb Arb::norm_pdf() const {
    return exp((this->pow(2)) / (-2)) / sqrt(2*pi());
}

Arb Arb::norm_pdf(const Arb& x) {
    return x.norm_pdf();
}

Arb Arb::norm_cdf() const {
    return (1 + ((*this)/sqrt(2)).erf())/2;
}

Arb Arb::norm_cdf(const Arb& x) {
    return x.norm_cdf();
}

Arb Arb::norm_cdf_inv() const {
    return sqrt(2)*erf_inv(2*(*this)-1);
}

Arb Arb::norm_cdf_inv(const Arb& x) {
    return x.norm_cdf_inv();
}


Arb Arb::sin() const {
    Arb ans;
    arb_sin(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}


Arb Arb::sin(const Arb& x) {
    return x.sin();
}


Arb Arb::acos() const {
    Arb ans;
    arb_acos(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}


Arb Arb::acos(const Arb& x) {
    return x.acos();
}