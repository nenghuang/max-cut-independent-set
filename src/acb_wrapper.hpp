#ifndef ACB_WRAPPER_HPP
#define ACB_WRAPPER_HPP

#include "flint/acb.h"
#include "arb_wrapper.hpp"

class Acb {
public:
    // constructors
    Acb();

    Acb(int real);

    Acb(double real);
    
    Acb(const Arb real);
    
    Acb(double real, double imag);

    Acb(const Arb real, const Arb imag);

    Acb(const acb_t z);
    
    // destructor
    virtual ~Acb();
    
    // debug/display info
    void println();

    // internal operators
    Acb operator+(const Acb& rhs) const;
    Acb operator-(const Acb& rhs) const;
    Acb operator*(const Acb& rhs) const;
    Acb operator/(const Acb& rhs) const;
  
    // interval arithmetic operations
    Arb real() const;
    Arb imag() const;

    int is_real() const;

    Acb join(const Acb& rhs) const;
    static Acb join(const Acb& lhs, const Acb& rhs);

    static Acb nan();
    int is_nan() const;

    // mathematical constants and functions
    static Acb pi();

    Acb exp() const;
    static Acb exp(const Acb& x);

    Acb sqrt() const;
    static Acb sqrt(const Acb& x);

    //useful for integration
    Acb sqrt_analytic(int analytic) const;
    static Acb sqrt_analytic(const Acb& x, int analytic); 

    // superior to pow(2)
    Acb sqr() const;
    static Acb sqr(const Acb& x);

    Acb pow(const Acb& rhs) const;
    static Acb pow(const Acb& lhs, const Acb& rhs);

    //useful for integration
    Acb pow_analytic(const Acb& rhs, int analytic) const;
    static Acb pow_analytic(const Acb& lhs, const Acb& rhs, int analytic);

    Acb erf() const;
    static Acb erf(const Acb& x);

    // following functions are for the standard Gaussian
    // with mean 0 and variance 1
    Acb norm_pdf() const;
    static Acb norm_pdf(const Acb& x);

    Acb norm_cdf() const;
    static Acb norm_cdf(const Acb& x);

    // internal data
    acb_t t;
};

// External operators
// Operations with a double on the LHS
// cannot be part of the class definition
Acb operator+(double lhs, const Acb& rhs);
Acb operator-(double lhs, const Acb& rhs);
Acb operator*(double lhs, const Acb& rhs);
Acb operator/(double lhs, const Acb& rhs);

// same with Arb on lhs
Acb operator+(const Arb& lhs, const Acb& rhs);
Acb operator-(const Arb& lhs, const Acb& rhs);
Acb operator*(const Arb& lhs, const Acb& rhs);
Acb operator/(const Arb& lhs, const Acb& rhs);

#endif
