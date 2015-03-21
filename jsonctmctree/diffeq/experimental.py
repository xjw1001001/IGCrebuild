"""
From the paper "Computing the Action of the Matrix Exponential":

If we wish to exponentiate the matrix t*A for several values of t then,
since alpha_p(t*A) = abs(t)*alph_p(A), we can precompute the matrix S
with pmax - 1 rows and mmax columns given by
S_p_m = alpha_p(A) / theta_m when
2 <= p <= pmax and p*(p-1)-1 <= m <= mmax,
and 0 otherwise,
and then for each t obtain C_mstar(t*A) as the smallest nonzero element
in the matrix ceil(abs(t)*S)*diag(1, 2, ..., mmax), where mstar
is the column index of the smallest element.
Table 3.1 lists somf of the theta_m values
corresponding to
u_s = tol = 2**-24 approx_eq 6e-8 (single precision) and
u_d = tol = 2**-53 approx_eq 1.1e-16 (double precision).
These values were determined as described in [14, App.].
That reference is Higham and Al-Mohy, Computing matrix functions (2010).

"""

