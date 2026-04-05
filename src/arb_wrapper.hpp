#ifndef ARB_WRAPPER_HPP
#define ARB_WRAPPER_HPP

#define GLOBAL_PRECISION 128 // Move this to a global file

#include "flint/arb.h"

class Arb {
public:
    // constructors
    Arb();
    
    Arb(double d); // 0-length interval around d
    
    Arb(double d1, double d2); // interval between d1 and d2
    
    // destructor
    virtual ~Arb();
    
    // debug/display info
    void print() const;
    void println() const;

    void pretty_print() const;
    void pretty_println() const;

    // export
    double left_edge_to_double() const;
    double mid_to_double() const;
    double right_edge_to_double() const;

    // internal operators
    Arb operator-() const; //unary negation

    Arb operator+(const Arb& rhs) const;
    Arb operator-(const Arb& rhs) const;
    Arb operator*(const Arb& rhs) const;
    Arb operator/(const Arb& rhs) const;

    Arb operator+(double rhs) const;
    Arb operator-(double rhs) const;
    Arb operator*(double rhs) const;
    Arb operator/(double rhs) const;

    // these interval operators are "for all"
    // x in lhs and y in rhs is op true?
    // in particular, == is meaningless unless
    // the radius of both are 0
    // also this means that == and != are *not*
    // complementary!
    int operator==(const Arb& rhs) const;
    int operator!=(const Arb& rhs) const;
    int operator<(const Arb& rhs) const;
    int operator<=(const Arb& rhs) const;
    int operator>(const Arb& rhs) const;
    int operator>=(const Arb& rhs) const;
  
    // interval arithmetic operations
    Arb join(const Arb& rhs) const;
    static Arb join(const Arb& lhs, const Arb& rhs);

    Arb intersect(const Arb& rhs) const;
    static Arb intersect(const Arb& lhs, const Arb& rhs);

    // checks if rhs inside lhs -- Note NaN treated like [-inf, +inf]
    int contains(const Arb& rhs) const;
    static int contains(const Arb& lhs, const Arb& rhs);

    // checks if lhs inside rhs
    int inside(const Arb& rhs) const;
    static int inside(const Arb& lhs, const Arb& rhs);

    Arb mid() const;
    Arb rad() const;

    Arb left_edge() const;
    Arb right_edge() const;

    Arb left_half() const;
    Arb right_half() const;

    static Arb nan();
    int is_nan() const;

    // mathematical constants and functions
    static Arb pi();

    Arb abs() const;
    static Arb abs(const Arb& x);

    Arb min(const Arb& rhs) const;
    static Arb min(const Arb& lhs, const Arb& rhs);

    Arb max(const Arb& rhs) const;
    static Arb max(const Arb& lhs, const Arb& rhs);

    Arb exp() const;
    static Arb exp(const Arb& x);

    Arb sqrt() const;
    static Arb sqrt(const Arb& x);

    // sqrt(max(x,0))
    Arb safe_sqrt() const;
    static Arb safe_sqrt(const Arb& x);

    Arb sqr() const;
    static Arb sqr(const Arb& x);

    Arb pow(const Arb& rhs) const;
    static Arb pow(const Arb& lhs, const Arb& rhs);

    Arb erf() const;
    static Arb erf(const Arb& x);

    Arb erf_inv() const;
    static Arb erf_inv(const Arb& x);

    // following functions are for the standard Gaussian
    // with mean 0 and variance 1
    Arb norm_pdf() const;
    static Arb norm_pdf(const Arb& x);

    Arb norm_cdf() const;
    static Arb norm_cdf(const Arb& x);

    Arb norm_cdf_inv() const;
    static Arb norm_cdf_inv(const Arb& x);

    Arb sin() const;
    static Arb sin(const Arb& rhs);

    Arb acos() const;
    static Arb acos(const Arb& rhs);

    // internal data -- should be private, but breaks acb_wrapper
    arb_t t;
};

// External operators
// Operations with a double on the LHS
// cannot be part of the class definition
Arb operator+(double lhs, const Arb& rhs);
Arb operator-(double lhs, const Arb& rhs);
Arb operator*(double lhs, const Arb& rhs);
Arb operator/(double lhs, const Arb& rhs);

#endif