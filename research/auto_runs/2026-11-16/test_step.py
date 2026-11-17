import pytest
import sympy as sp
import numpy as np


def A_m_val(m_val, gamma_val):
    """Evaluate A_m(gamma) = m(m+3)/(2*gamma) - m(m+1)/gamma**2."""
    return m_val * (m_val + 3) / (2 * gamma_val) - m_val * (m_val + 1) / gamma_val**2


def partial_sum_val(M_val, gamma_val):
    """Evaluate closed-form candidate for sum_{m=0}^{M} A_m(gamma)."""
    return M_val * (M_val + 1) / 3 * (
        (M_val + 5) / (2 * gamma_val) - (M_val + 2) / gamma_val**2
    )


class TestAmKnownIdentities:
    """Re-verify the two Stage-60 identities as baseline."""

    def test_A0_vanishes(self):
        m, g = sp.symbols('m gamma', positive=True)
        A = m * (m + 3) / (2 * g) - m * (m + 1) / g**2
        assert sp.simplify(A.subs(m, 0)) == 0

    def test_minus_gamma_A1_equals_c_gamma(self):
        g = sp.Symbol('gamma', positive=True)
        A1 = 1 * 4 / (2 * g) - 1 * 2 / g**2
        c_gamma = 2 * (1 - g) / g
        assert sp.simplify(-g * A1 - c_gamma) == 0


class TestAmSignDefiniteness:
    """A_m(gamma) < 0 for all m >= 1, gamma in (0,1)."""

    def test_sign_threshold_symbolic(self):
        """A_m < 0 iff gamma < 2(m+1)/(m+3); verify threshold >= 1 for m >= 1."""
        m = sp.Symbol('m', positive=True)
        threshold = 2 * (m + 1) / (m + 3)
        # threshold - 1 = (m-1)/(m+3) >= 0 for m >= 1
        gap = sp.simplify(threshold - 1)
        assert sp.simplify(gap - (m - 1) / (m + 3)) == 0
        # At m=1, gap=0 so threshold=1 exactly; gamma < 1 strictly => A_1 < 0
        assert gap.subs(m, 1) == 0
        assert gap.subs(m, 2) > 0

    def test_sign_numeric_grid(self):
        """Brute-force: A_m(gamma) < 0 on a 100x100 grid."""
        gamma_vals = np.linspace(0.01, 0.999, 100)
        m_vals = np.arange(1, 101)
        for gv in gamma_vals:
            vals = m_vals * (m_vals + 3) / (2 * gv) - m_vals * (m_vals + 1) / gv**2
            assert np.all(vals < -1e-15), (
                f"Non-negative A_m found at gamma={gv:.4f}"
            )


class TestPartialSumClosedForm:
    """Verify closed form of sum_{m=0}^{M} A_m(gamma)."""

    def test_closed_form_symbolic(self):
        m, g, M = sp.symbols('m gamma M', positive=True)
        A_expr = m * (m + 3) / (2 * g) - m * (m + 1) / g**2
        computed = sp.summation(A_expr, (m, 0, M))
        expected = M * (M + 1) / 3 * ((M + 5) / (2 * g) - (M + 2) / g**2)
        assert sp.simplify(sp.expand(computed) - sp.expand(expected)) == 0

    def test_closed_form_numeric(self):
        """Compare closed form vs direct summation."""
        for gv in [0.1, 0.3, 0.5, 0.7, 0.9]:
            for Mi in [1, 5, 10, 25, 50, 100]:
                direct = sum(A_m_val(m, gv) for m in range(Mi + 1))
                closed = partial_sum_val(Mi, gv)
                assert abs(direct - closed) < 1e-10 * (abs(direct) + 1), (
                    f"Mismatch at M={Mi}, gamma={gv}: "
                    f"direct={direct}, closed={closed}"
                )

    def test_partial_sum_sign(self):
        """Partial sum is negative for all M >= 1, gamma in (0,1)."""
        for Mi in range(1, 61):
            for gv in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.999]:
                val = partial_sum_val(Mi, gv)
                assert val < 0, f"Positive partial sum at M={Mi}, gamma={gv}: {val}"


class TestInnerLayerScaling:
    """The O(1/n) inner-layer contribution scales as Theta(M^3/n)."""

    def test_leading_coefficient_negative(self):
        """Leading coeff of partial sum in M is (gamma-2)/(6*gamma^2) < 0 for gamma in (0,1)."""
        g = sp.Symbol('gamma', positive=True)
        leading = (g - 2) / (6 * g**2)
        # Verify this matches the M^3 coefficient of the closed form
        M = sp.Symbol('M', positive=True)
        closed = M * (M + 1) / 3 * ((M + 5) / (2 * g) - (M + 2) / g**2)
        poly = sp.Poly(sp.expand(closed), M)
        # Extract M^3 coefficient
        coeffs = poly.all_coeffs()  # highest degree first
        m3_coeff = sp.simplify(coeffs[0])  # coefficient of M^3
        assert sp.simplify(m3_coeff - leading) == 0
        # Verify negative for gamma in (0,2)
        for gv in [0.1, 0.5, 0.9, 1.5]:
            assert float(leading.subs(g, gv)) < 0

    def test_scaling_hierarchy(self):
        """
        Verify contribution = |partial_sum / n| obeys:
        - M = n^{1/4}: contribution -> 0 as n -> infty
        - M = n^{1/2}: contribution -> infty as n -> infty
        
        This proves the inner expansion alone cannot supply
        the O(1) crossover mass; it must come from m ~ sqrt(n).
        """
        for gv in [0.3, 0.5, 0.8]:
            contribs_quarter = []
            contribs_half = []
            for exp_n in [4, 6, 8, 10]:
                n = 10**exp_n
                Mq = max(2, int(round(n**0.25)))
                Mh = max(2, int(round(n**0.5)))
                cq = abs(partial_sum_val(Mq, gv) / n)
                ch = abs(partial_sum_val(Mh, gv) / n)
                contribs_quarter.append(cq)
                contribs_half.append(ch)

            # Quarter-power contributions should DECREASE
            for i in range(1, len(contribs_quarter)):
                assert contribs_quarter[i] < contribs_quarter[i - 1], (
                    f"gamma={gv}: n^(1/4) contribs not decreasing: {contribs_quarter}"
                )
            # Half-power contributions should INCREASE
            for i in range(1, len(contribs_half)):
                assert contribs_half[i] > contribs_half[i - 1], (
                    f"gamma={gv}: n^(1/2) contribs not increasing: {contribs_half}"
                )
            # At n=10^10, quarter-power should be small
            assert contribs_quarter[-1] < 1.0, (
                f"gamma={gv}: n^(1/4) contrib at n=10^10 too large: {contribs_quarter[-1]}"
            )
            # At n=10^10, half-power should be large
            assert contribs_half[-1] > 100, (
                f"gamma={gv}: n^(1/2) contrib at n=10^10 too small: {contribs_half[-1]}"
            )
