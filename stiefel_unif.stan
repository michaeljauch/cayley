functions {
  // Compute the Kronecker product. 
    matrix kron(matrix A, matrix B) {
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron;
    for (i in 1:cols(A)) {
      for (j in 1:rows(A)) {
        kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = 
        A[j,i] * B;
      }
    }
    return kron;
  }
  // Calculate the log of the m-dimensional Jacobian Jm according to the 
  // approach described in the appendix of the paper. No terms are dropped due 
  // to proportionality.
  real calc_logJm(int p, int k, vector phi, matrix Dk, matrix K) {
  int m; 
  matrix[k,k] B; 
  matrix[p-k,k] A;
  matrix[k,k] C11; 
  matrix[k, p-k] C12; 
  matrix[p-k, k] C21; 
  matrix[p-k, p-k] C22;
  matrix[k,k] G11; 
  matrix[k, p-k] G12; 
  matrix[p-k, k] G21; 
  matrix[p-k, p-k] G22;
  matrix[k,k] H11; 
  matrix[k, p-k] H12; 
  matrix[p-k, k] H21; 
  matrix[p-k, p-k] H22;
  matrix[k*(k-1)/2, k*(k-1)/2] Om11;
  matrix[k*(p-k), k*(p-k)] Om22; 
  matrix[k*(k-1)/2, k*(p-k)] Om12; 
  matrix[k,k] Ik; 
  matrix[p-k, p-k] Ipmink; 
  matrix[k*(p-k), k*(p-k)] L22; // Cholesky factor of Om22
  real logJm;
  m = p*k -k*(k+1)/2;
  Ik = diag_matrix(rep_vector(1.0,k)); 
  Ipmink = diag_matrix(rep_vector(1.0,p-k)); 
  B = to_matrix(Dk*phi[1:k*(k-1)/2], k, k); 
  A = to_matrix(phi[k*(k-1)/2+1:p*k - k*(k+1)/2], p-k, k); 
  C11 = inverse(Ik - B + A'*A); 
  C12 = -C11*A'; 
  C21 = A*C11; 
  C22 = Ipmink - A*C11*A';
  H11 = C11'*C11 + C21'*C21; 
  H12 = C11'*C12 + C21'*C22; 
  H21 = C12'*C11 + C22'*C21; 
  H22 = C12'*C12 + C22'*C22; 
  G11 = C11*C11'; 
  G12 = C11*C21'; 
  G21 = C21*C11'; 
  G22 = C21*C21'; 
  Om11 = Dk'*kron(G11, H11)*Dk;
  Om22 = kron(G11, H22) + kron(H11, G22) - (kron(G12, H21) + kron(H12, G21))*K;
  L22 = cholesky_decompose(Om22);
  Om12 = Dk'*kron(G11, H12) - Dk'*kron(G12, H11)*K;
  logJm = m*log(2) + .5*(2*sum(log(diagonal(L22)))) + 
  .5*log_determinant(Om11 - crossprod(mdivide_left_tri_low(L22, Om12')));
  return logJm;
  }
                                                           
}

data{
int<lower=1> k; // of columns of Stiefel manifold elt.
int<lower=k+1> p; // of rows of Stiefel manifold elt. 
matrix[k*k, k*(k-1)/2] Dk; // D*_k matrix
matrix[k*(p-k), k*(p-k)] K; // Commutation matrix K_p-k,k
}

parameters{
  vector[p*k - k*(k+1)/2] phi;
}

transformed parameters{
matrix[p,k] Q; 
{
  matrix[k,k] B; 
  matrix[p-k,k] A;
  matrix[k,k] F;
  matrix[p, k] LHS; 
  matrix[k, k] RHS;
  B = to_matrix(Dk*phi[1:k*(k-1)/2], k, k); 
  A = to_matrix(phi[k*(k-1)/2+1:p*k - k*(k+1)/2], p-k, k); 
  F = A'*A - B; 
  LHS = append_row(diag_matrix(rep_vector(1.0,k)) - F, 2*A);
  RHS = inverse(diag_matrix(rep_vector(1.0,k)) + F); 
  Q = LHS*RHS; 
}
}

model{
  target += calc_logJm(p, k, phi, Dk, K);
}
