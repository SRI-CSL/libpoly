set(poly_SOURCES 
    polynomial/polynomial.c
    polynomial/polynomial_internal.c
    polynomial/polynomial_context.c
    polynomial/coefficient.c
    polynomial/monomial.c
    upolynomial/bounds.c
    upolynomial/factorization.c
    upolynomial/upolynomial.c
    upolynomial/gcd.c
    upolynomial/upolynomial_dense.c
    upolynomial/root_finding.c
    upolynomial/upolynomial_internal.c
    number/rational.c
    number/algebraic_number.c
    number/dyadic_rational.c
    number/integer.c
    variable/variable.c
    variable/variable_order.c
    interval/interval.c
    utils/debug_trace.c
    utils/assignment_internal.c
    utils/assignment.c
    utils/statistics.c
)

add_library(poly ${poly_SOURCES})
