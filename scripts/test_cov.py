#!/usr/bin/env python3

# Debug covariance calculations using NumPy

import numpy as np
import scipy

w_file = "DEBUG_w.out"
n_file = "DEBUG_N.out"

w     = np.loadtxt(w_file, comments='#', dtype=np.float64)
n_vec = np.loadtxt(n_file, comments='#', dtype=np.float64)

num_samples_total = w.shape[0]  # N
num_states        = w.shape[1]  # K

n    = np.diag(n_vec)
wT_w = np.matrix(w.T @ w)

# Check properties of 'w'
# - Rows
sum_over_rows = np.sum(w, axis=0)
err_rowsum = np.linalg.norm(sum_over_rows - 1.0)
# - Columns
sum_over_columns = np.matrix( np.empty((num_samples_total,1), dtype=np.float64) )
for r in range(0, num_samples_total):
	sum_over_columns[r] = np.inner( n_vec, w[r,:] )
err_colsum = np.linalg.norm(sum_over_columns - 1.0)
print("Check the matrix of weights")
print("  colsum:  sum_{n=1}^N w_{n,i} = 1       -->  total_error = ", err_rowsum)
print("  rowsum:  sum_{i=1}^K N_i * w_{n,i} = 1 -->  total_error = ", err_colsum)

# SVD
def svd():
	u, s, vT = np.linalg.svd(w, full_matrices=True, compute_uv=True, hermitian=False)
	s[np.where(s < 0.0)] = 0.0
	sigma = np.matrix( np.diag(s) )
	v     = np.matrix( vT.T )
	return sigma, v

# Eigenvalue decomposition
def eig():
	lambdas, v = np.linalg.eigh(wT_w)
	lambdas[np.where(lambdas < 0.0)] = 0.0
	sigma = np.matrix( np.diag(np.sqrt(lambdas)) )
	v     = np.matrix( v )
	return sigma, v


# Run the chosen method
#method = svd
method = eig
sigma, v = method()
vT = np.matrix(v.T)
#vT = np.matrix( np.transpose(v) )

print("\nSingular values: sigma =\n", sigma)
print("type(w)     =", type(w))
print("type(wT_w)  =", type(wT_w))
print("type(sigma) =", type(sigma))
print("type(v)     =", type(v))

# Pseudo-inverse
v_sigma = np.matrix( v @ sigma )
i_K     = np.matrix( np.eye(num_states) )
p       = v_sigma.T @ n @ v_sigma
m       = np.matrix( i_K - p )
#m       = np.matrix( i_K - (sigma @ v.T) @ n @ (v @ sigma) )
pinv_m  = np.matrix( np.linalg.pinv(m) )
print("\np = \n", p)
print("  cond(p)   =", np.linalg.cond(p))
print("\nm = \n", m)
print("\nm^+ = ", pinv_m)
print("  rank(m)   =", np.linalg.matrix_rank(m))
print("  cond(m)   =", np.linalg.cond(m))
print("  cond(m^+) =", np.linalg.cond(pinv_m))
print("  err(m^+)  = ||m*(m^+)*m - m||_F =", np.linalg.norm(m*pinv_m*m - m, ord='fro'))
print("\nv*sigma = ", v_sigma)
print("  cond(v*sigma) =", np.linalg.cond(v_sigma))

# Covariance matrix and corresponding standard error (from diagonal)
max_abs = np.max(np.abs(pinv_m))
cov_f = v_sigma @ pinv_m @ v_sigma.T
#cov_f = v_sigma @ pinv_m @ v_sigma.T
#cov_f = (v @ sigma) @ pinv_m @ (sigma @ v.T)
std_err_f = np.sqrt( np.diagonal(cov_f) )
print("\ncov_f = ")
print(cov_f)
print("std_err_f = ")
print(std_err_f)

# TEST: Approximation from Kong et al. 
print("\nALT: Approximate Theta ~ W^T*W")
cov_f_alt = w.T @ w
std_err_f_alt = np.sqrt( np.diagonal(cov_f_alt) )
print("std_err_f(alt) =")
print(std_err_f_alt)
print("")

# TEST: Pseudo-inverse via SVD of 'm'
print("\nALT: Directly compute M^+ via SVD")
u_m, s_m, vT_m = np.linalg.svd(m, full_matrices=True, compute_uv=True, hermitian=False)
s_m[np.where(s_m < 0.0)] = 0.0
sigma_m = np.matrix( np.diag(s_m) )
v_m     = np.matrix( vT_m.T )
# - Invert Sigma_M
tol = 1.0e-10
inv_sigma_m = np.matrix( np.zeros((num_states, num_states)), dtype=np.float64)
for i in range(0, num_states):
	if ( sigma_m[i,i] > tol ):
		inv_sigma_m[i,i] = 1.0/sigma_m[i,i]
# - Invert M
pinv_m_alt = v_m @ inv_sigma_m @ u_m.T
# - Finish up
cov_f_alt  = v_sigma @ pinv_m_alt @ v_sigma.T
std_err_f_alt = np.sqrt( np.diagonal(cov_f_alt) )
print("pinv_m(alt) =")
print(pinv_m_alt)
print("cov_f(alt) =")
print(cov_f_alt)
print("std_err_f(alt) =")
print(std_err_f_alt)
print("err(diff. M^+) =", np.linalg.norm(np.abs(pinv_m_alt - pinv_m)))


### Scraps ###

# This uses an incredible amount of memory and takes forever
#print("\nALT: Use W directly")
#m_alt = np.eye(num_samples_total) - (w @ n @ w.T)
#cov_f_alt = w.T @ np.linalg.pinv(m_alt) @ w
#std_err_f_alt = np.sqrt( np.diagonal(cov_f_alt) )
#print("std_err_f(alt) = ")
#print(std_err_f_alt)

# Comparing methods would require demuxing the output...
#sigma_svd, v_svd = svd()
#sigma_eig, v_eig = eig()
#print("||Delta v||_F     =", np.linalg.norm(v_svd     - v_eig,     ord='fro'))
#print("||Delta sigma||_F =", np.linalg.norm(sigma_svd - sigma_eig, ord='fro'))
