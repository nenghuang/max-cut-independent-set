#include "acb_wrapper.hpp"
#include "flint/acb_hypgeom.h"
#include <cassert>

// constructors

Acb::Acb() {
    acb_init(this -> t);
}

Acb::Acb(int real) {
    acb_init(this -> t);
    acb_set_si(this -> t, real);
}


Acb::Acb(double real) {
    acb_init(this -> t);
    acb_set_d(this -> t, real);
}

Acb::Acb(const Arb real) {
    acb_init(this -> t);
    acb_set_arb(this -> t, real.t);
}

Acb::Acb(double real, double imag) {
    acb_init(this -> t);
    acb_set_d_d(this -> t, real, imag);
}

Acb::Acb(const Arb real, const Arb imag) {
    acb_init(this -> t);
    acb_set_arb_arb(this -> t, real.t, imag.t);
}

Acb::Acb(const acb_t z) {
    acb_init(this -> t);
    acb_set(this -> t, z);
}

// destructor

Acb::~Acb() {
    acb_clear(this -> t);
}

// display info

void Acb::println() {
    acb_printd(this -> t, 10);
    flint_printf("\n");
}

// interval arithmetic operations
Arb Acb::real() const {
    Arb ans;
    acb_get_real(ans.t, this -> t);
    return ans;
}

Arb Acb::imag() const {
    Arb ans;
    acb_get_imag(ans.t, this -> t);
    return ans;
}

int Acb::is_real() const {
    return acb_is_real(this -> t);
}


Acb Acb::join(const Acb& rhs) const {
    Acb ans;
    acb_union(ans.t, this -> t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::join(const Acb& lhs, const Acb& rhs) {
    return lhs.join(rhs);
}

Acb Acb::nan() {
    Acb ans;
    acb_indeterminate(ans.t);
    return ans;
}

int Acb::is_nan() const {
    return (this->real().is_nan()) || (this->imag().is_nan());
}

// operators

Acb Acb::operator+(const Acb& rhs) const {
    Acb ans;
    acb_add(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::operator-(const Acb& rhs) const {
    Acb ans;
    acb_sub(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::operator*(const Acb& rhs) const {
    Acb ans;
    acb_mul(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::operator/(const Acb& rhs) const {
    Acb ans;
    acb_div(ans.t, this->t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Acb operator+(double lhs, const Acb& rhs) {
    Acb x(lhs);
    return x + rhs;
}

Acb operator-(double lhs, const Acb& rhs) {
    Acb x(lhs);
    return x - rhs;
}

Acb operator*(double lhs, const Acb& rhs) {
    Acb x(lhs);
    return x * rhs;
}

Acb operator/(double lhs, const Acb& rhs) {
    Acb x(lhs);
    return x / rhs;
}

Acb operator+(const Arb& lhs, const Acb& rhs) {
    Acb x(lhs);
    return x + rhs;
}

Acb operator-(const Arb& lhs, const Acb& rhs) {
    Acb x(lhs);
    return x - rhs;
}

Acb operator*(const Arb& lhs, const Acb& rhs) {
    Acb x(lhs);
    return x * rhs;
}

Acb operator/(const Arb& lhs, const Acb& rhs) {
    Acb x(lhs);
    return x / rhs;
}

// mathematical constants and functions

Acb Acb::pi() {
    Acb ans;
    acb_const_pi(ans.t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::exp() const {
    Acb ans;
    acb_exp(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::exp(const Acb& x) {
    return x.exp();
}

Acb Acb::sqrt() const {
    Acb ans;
    acb_sqrt(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::sqrt(const Acb& x) {
    return x.sqrt();
}

Acb Acb::sqrt_analytic(int analytic) const {
    Acb ans;
    acb_sqrt_analytic(ans.t, this -> t, analytic, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::sqrt_analytic(const Acb& x, int analytic) {
    return x.sqrt_analytic(analytic);
}

Acb Acb::sqr() const {
    Acb ans;
    acb_sqr(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::sqr(const Acb& x) {
    return x.sqr();
}

Acb Acb::pow(const Acb& rhs) const {
    Acb ans;
    acb_pow(ans.t, this -> t, rhs.t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::pow(const Acb& lhs, const Acb& rhs) {
    return lhs.pow(rhs);
}

Acb Acb::pow_analytic(const Acb& rhs, int analytic) const {
    Acb ans;
    acb_pow_analytic(ans.t, this -> t, rhs.t, analytic, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::pow_analytic(const Acb& lhs, const Acb& rhs, int analytic) {
    return lhs.pow_analytic(rhs, analytic);
}

Acb Acb::erf() const {
    Acb ans;
    acb_hypgeom_erf(ans.t, this -> t, GLOBAL_PRECISION);
    return ans;
}

Acb Acb::erf(const Acb& x) {
    return x.erf();
}

Acb Acb::norm_pdf() const {
    return exp((this->sqr()) / (-2)) / sqrt(2*pi());
}

Acb Acb::norm_pdf(const Acb& x) {
    return x.norm_pdf();
}

Acb Acb::norm_cdf() const {
    return (1 + ((*this)/sqrt(2)).erf())/2;
}

Acb Acb::norm_cdf(const Acb& x) {
    return x.norm_cdf();
}
