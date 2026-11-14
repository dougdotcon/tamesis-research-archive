import pytest
import numpy as np
from sympy import symbols, Rational, simplify, series, S


gamma = symbols('gamma', positive=True)
m_s = symbols('m', positive=True, integer=True)


def A(m_val, g):
    """A_m(gamma) = m(m+3)/(2*gamma) - m(m+1)/gamma**2."""
    return m_val * (m_val + 3) / (2 * g) - m_val * (m_val + 1) / g**2


# ── 1. Established identities ──────────────────────────────────────

class TestEstablishedIdentities:
    def test_A0_identically_zero(self):
        """A_0(gamma) = 0 for all gamma."""
        assert A(0, gamma) == 0

    def test_minus_gamma_A1_equals_c_gamma(self):
        """-gamma * A_1(gamma) = c(gamma) = 2(1-gamma)/gamma."""
        lhs = simplify(-gamma * A(1, gamma))
        rhs = 2 * (1 - gamma) / gamma
        assert simplify(lhs - rhs) == 0


# ── 2. Sign structure ──────────────────────────────────────────────

class TestSignStructure:
    def test_gamma_crit_geq_1_symbolic(self):
        """gamma_crit(m) = 2(m+1)/(m+3) >= 1 for m >= 1.
        Proof: 2(m+1) - (m+3) = m - 1 >= 0 for m >= 1."""
        expr = simplify(2 * (m_s + 1) - (m_s + 3))
        assert simplify(expr - (m_s - 1)) == 0

    @pytest.mark.parametrize("mm", range(1, 25))
    def test_Am_negative_on_unit_interval(self, mm):
        """A_m(gamma) < 0 for integer m >= 1, gamma in (0,1)."""
        for g in np.linspace(0.05, 0.95, 30):
            val = mm * (mm + 3) / (2.0 * g) - mm * (mm + 1) / g**2
            assert val < 0, f"A_{mm}({g:.2f}) = {val:.6f} >= 0"


# ── 3. NEW: Partial-sum closed form ────────────────────────────────

class TestPartialSumClosedForm:
    """Derives and verifies:
    sum_{m=0}^{M} A_m(gamma) = M(M+1)[gamma(M+5) - 2(M+2)] / (6 gamma^2).
    """

    @pytest.mark.parametrize("M_val", range(0, 20))
    def test_partial_sum_vs_direct(self, M_val):
        direct = sum(A(i, gamma) for i in range(M_val + 1))
        formula = (S(M_val) * (M_val + 1)
                   * (gamma * (M_val + 5) - 2 * (M_val + 2))
                   / (6 * gamma**2))
        assert simplify(direct - formula) == 0, f"Mismatch at M={M_val}"

    def test_partial_sum_negative_for_M_geq_1(self):
        """P_M(gamma) < 0 for M >= 1, gamma in (0,1):
        bracket gamma*(M+5) - 2*(M+2) = 0 at gamma=2(M+2)/(M+5) >= 1
        for M >= 1."""
        for M_val in range(1, 15):
            gc = 2.0 * (M_val + 2) / (M_val + 5)
            assert gc >= 1.0, f"gamma_crit({M_val}) = {gc} < 1"


# ── 4. NEW: Generating function ────────────────────────────────────

class TestGeneratingFunction:
    """G(x, gamma) = x * [(2-x)*gamma - 2] / [gamma^2 * (1-x)^3]
    satisfies: coeff(x^m, G) = A_m(gamma) for all m >= 0.
    """

    def test_gf_coefficients_match_Am(self):
        x = symbols('x')
        gf = x * ((2 - x) * gamma - 2) / (gamma**2 * (1 - x)**3)
        s = series(gf, x, 0, n=13)
        for mm in range(13):
            coeff = s.coeff(x, mm)
            expected = A(mm, gamma)
            assert simplify(coeff - expected) == 0, (
                f"GF coeff at m={mm}: got {coeff}, expected {expected}")


# ── 5. Second difference constant in m ─────────────────────────────

class TestSecondDifference:
    def test_delta2_Am_independent_of_m(self):
        """Delta^2 A_m = (gamma - 2) / gamma^2, independent of m.
        (Consequence of A_m being exactly quadratic in m.)"""
        d2 = A(m_s + 2, gamma) - 2 * A(m_s + 1, gamma) + A(m_s, gamma)
        assert simplify(d2 - (gamma - 2) / gamma**2) == 0


# ── 6. Exact rational spot-checks ──────────────────────────────────

class TestExactRationalValues:
    @pytest.mark.parametrize("mm,g,expected", [
        (0, Rational(1, 2), 0),
        (1, Rational(1, 2), -4),
        (2, Rational(1, 2), -14),
        (3, Rational(1, 2), -30),
        (1, Rational(1, 3), -12),
        (2, Rational(1, 3), -39),
        (3, Rational(1, 3), -90),
        (1, Rational(4, 5), Rational(-5, 8)),
    ])
    def test_exact_value(self, mm, g, expected):
        assert A(mm, g) == expected, (
            f"A_{mm}({g}) = {A(mm, g)}, expected {expected}")
