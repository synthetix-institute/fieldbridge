"""Finite-count inference for paired changes in occupation of a stored state."""
from __future__ import annotations

import math

from scipy.stats import beta, binomtest


def contrast(plus, minus, n, alpha=.05):
    """Exact conditional sign test and conservative marginal-binomial interval.

    The two discordant proportions form a multinomial, so four one-sided
    Clopper-Pearson bounds use alpha/4 each. Their difference has simultaneous
    coverage of at least 1-alpha, including when no changes are observed.
    """
    if any(int(k)!=k or k<0 for k in (plus,minus,n)) or n<2 or plus+minus>n or not 0<alpha<1:
        raise ValueError("Invalid paired occupation counts")
    plus,minus,n=map(int,(plus,minus,n))
    tail=alpha/4
    def interval(k):
        return (float(beta.ppf(tail,k,n-k+1)) if k else 0.,
                float(beta.ppf(1-tail,k+1,n-k)) if k<n else 1.)
    lp,up=interval(plus)
    lm,um=interval(minus)
    effect=(plus-minus)/n
    variance=max(0.,((plus+minus)-n*effect**2)/(n-1))
    return {"discordant_plus":plus,"discordant_minus":minus,"n":n,"effect":effect,
            "standard_error":math.sqrt(variance/n),"lower_95":lp-um,"upper_95":up-lm,
            "p_value":float(binomtest(plus,plus+minus,.5).pvalue) if plus+minus else 1.,
            "interval":"simultaneous exact marginal-binomial bounds"}


def from_summary(effect,standard_error,n):
    """Recover the integer discordances from saved mean and paired sample SE."""
    difference=round(effect*n)
    total=round(n*(n-1)*standard_error**2+n*effect**2)
    if (total+difference)%2 or abs(effect*n-difference)>1e-6:
        raise ValueError("The summary does not describe paired binary outcomes")
    return contrast((total+difference)//2,(total-difference)//2,n)
