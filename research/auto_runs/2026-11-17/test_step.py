import pytest
import sympy as sp
import numpy as np


def test_A_m_polynomial_decomposition():
    """Verify A_m(g) = alpha(g)*m^2 + beta(g)*m."""
    m = sp.Symbol('m', nonnegative=True, integer=True)
    g = sp.Symbol('gamma', positive=True)
    A = m * (m + 3) / (2 * g) - m * (m + 1) / g**2
    alpha = (g - 2) / (2 * g**2)
    beta = (3 * g - 2) / (2 * g**2)
    assert sp.simplify(A - (alpha * m**2 + beta * m)) == 0


def test_beta_unique_zero():
    """beta(g) = (3g-2)/(2g^2) vanishes uniquely at g=2/3 in (0,1)."""
    g = sp.Symbol('gamma', positive=True)
    beta = (3 * g - 2) / (2 * g**2)
    zeros = sp.solve(beta, g)
    assert zeros == [sp.Rational(2, 3)]


def test_A_m_at_gamma_two_thirds():
    """At g=2/3, A_m = -3m^2/2 (pure quadratic, zero linear term)."""
    m = sp.Symbol('m', nonnegative=True, integer=True)
    g = sp.Symbol('gamma', positive=True)
    A = m * (m + 3) / (2 * g) - m * (m + 1) / g**2
    A_23 = A.subs(g, sp.Rational(2, 3))
    expected = sp.Rational(-3, 2) * m**2
    assert sp.simplify(A_23 - expected) == 0, f"got {A_23}"
    for mv in range(51):
        val = float(A_23.subs(m, mv))
        assert abs(val - (-1.5 * mv**2)) < 1e-14, f"m={mv}"


def test_c_gamma_consistency_at_two_thirds():
    """-g*A_1(g) = c(g) = 2(1-g)/g; at g=2/3 both equal 1."""
    g = sp.Symbol('gamma', positive=True)
    m = sp.Symbol('m', nonnegative=True, integer=True)
    A = m * (m + 3) / (2 * g) - m * (m + 1) / g**2
    c_g = 2 * (1 - g) / g
    assert sp.simplify(-g * A.subs(m, 1) - c_g) == 0
    assert float(c_g.subs(g, sp.Rational(2, 3))) == 1.0
    assert float((-g * A.subs(m, 1)).subs(g, sp.Rational(2, 3))) == 1.0


def test_strict_negativity():
    """A_m(g) < 0 for all m >= 1 and g in (0,1).

    Proof sketch: A_m = m/(2g^2)*[(m+3)g - 2(m+1)], zero at
    g* = 2(m+1)/(m+3). Since g*-1 = (m-1)/(m+3) >= 0 for m>=1,
    g* >= 1 > g, so the bracket is negative."""
    m = sp.Symbol('m', nonnegative=True, integer=True)
    g = sp.Symbol('gamma', positive=True)
    A = m * (m + 3) / (2 * g) - m * (m + 1) / g**2
    gstar = 2 * (m + 1) / (m + 3)
    assert sp.simplify(A.subs(g, gstar)) == 0
    assert sp.simplify(gstar - 1 - (m - 1) / (m + 3)) == 0
    for mv in range(1, 30):
        for gv in np.linspace(0.01, 0.99, 50):
            Av = mv * (mv + 3) / (2 * gv) - mv * (mv + 1) / gv**2
            assert Av < -1e-15, f"A_{mv}({gv:.3f}) = {Av}"


def test_second_difference_constant():
    """Delta^2 A_m = A_{m+2}-2A_{m+1}+A_m = (g-2)/g^2, constant in m."""
    m = sp.Symbol('m', nonnegative=True, integer=True)
    g = sp.Symbol('gamma', positive=True)
    def A(k):
        return k * (k + 3) / (2 * g) - k * (k + 1) / g**2
    delta2 = A(m + 2) - 2 * A(m + 1) + A(m)
    expected = (g - 2) / g**2
    assert sp.simplify(delta2 - expected) == 0
    assert float(expected.subs(g, sp.Rational(2, 3))) == -3.0


def test_monotonicity_transition():
    """A_m/m^2 = alpha+beta/m approaches alpha from below when g<2/3,
    from above when g>2/3, constant at -3/2 when g=2/3."""
    # g < 2/3: beta < 0 -> increasing toward limit from below
    for gv in [0.1, 0.3, 0.5, 0.6]:
        alpha_v = (gv - 2) / (2 * gv**2)
        ratios = []
        for mv in range(1, 40):
            Av = mv * (mv + 3) / (2 * gv) - mv * (mv + 1) / gv**2
            ratios.append(Av / mv**2)
        for i in range(len(ratios) - 1):
            assert ratios[i + 1] > ratios[i] - 1e-14, \
                f"g={gv}: not increasing at m={i+2}"
        for r in ratios:
            assert r < alpha_v + 1e-14, f"g={gv}: exceeds limit"
    # g > 2/3: beta > 0 -> decreasing toward limit from above
    for gv in [0.7, 0.8, 0.9, 0.95]:
        alpha_v = (gv - 2) / (2 * gv**2)
        ratios = []
        for mv in range(1, 40):
            Av = mv * (mv + 3) / (2 * gv) - mv * (mv + 1) / gv**2
            ratios.append(Av / mv**2)
        for i in range(len(ratios) - 1):
            assert ratios[i + 1] < ratios[i] + 1e-14, \
                f"g={gv}: not decreasing at m={i+2}"
        for r in ratios:
            assert r > alpha_v - 1e-14, f"g={gv}: below limit"
    # g = 2/3: constant
    for mv in range(1, 40):
        gv = 2.0 / 3.0
        Av = mv * (mv + 3) / (2 * gv) - mv * (mv + 1) / gv**2
        assert abs(Av / mv**2 - (-1.5)) < 1e-13, f"m={mv}"


def test_A0_vanishes():
    """A_0(g) = 0 identically (re-verified)."""
    g = sp.Symbol('gamma', positive=True)
    A0 = 0 * (0 + 3) / (2 * g) - 0 * (0 + 1) / g**2
    assert A0 == 0


def test_partial_sum_closed_form():
    """Sum_{m=0}^{M} A_m = alpha*M(M+1)(2M+1)/6 + beta*M(M+1)/2."""
    for gv in [0.3, 0.5, 2.0/3.0, 0.8]:
        alpha_v = (gv - 2) / (2 * gv**2)
        beta_v = (3 * gv - 2) / (2 * gv**2)
        for M in range(1, 30):
            direct = sum(
                mv * (mv + 3) / (2 * gv) - mv * (mv + 1) / gv**2
                for mv in range(M + 1)
            )
            formula = (alpha_v * M * (M + 1) * (2 * M + 1) / 6
                       + beta_v * M * (M + 1) / 2)
            assert abs(direct - formula) < 1e-9, \
                f"g={gv}, M={M}: {direct} vs {formula}"
