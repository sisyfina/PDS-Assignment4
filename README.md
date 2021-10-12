# Exercise 4 - Parallel & Distributed Computer Systems

Sparse Boolean Matrix Multiplication with two levels of parallelization

In this assignment we will implement Boolean matrix multiplication (BMM) of sparse CSC matrices using at least two forms of parallelism that combined produce the correct result faster than either form of parallelism independently. You may combine shared and distributed memory parallelism either OpenMP or Cilk, together with MPI, and or GPU programming.

We denote the logical conjunction (and), disjunction (or) and negation (not) by ∧, ∨ and ¬, respectively. Given matrix A = A(I, K) = [a(i, k)] with row index set I and column index set K, matrix B = B(K, J) = [b(k, j)] with index sets K and J. Let ni = |I|, nj = |J| and nk = |K|. Let C = C(I, J) = [c(i, j)] be the Boolean product C = F ⊙ (AB) c(i, j) = ⋁k∈K a(i, k) ∧ b(k, j), k ∈ K, i ∈ I, j ∈ J, F(i, j) ≠ 0 

Element c(i,j) = 1 if and only if there exists k in the search domain K such that a(i,k)b(k,j) = 1, i.e., there is a match between the binary vector in row i of A and the binary vector in column j of B.

If provided a filter F on the product AB, then only those elements with F(i,j) = 1 are calculated. The filtered product can be expressed as F ⊙ (AB), where ⊙ denotes the Hadamard product, i.e., element-wise product. The absence of a prescribed filter is equivalent to the case where F has no zero elements. We may therefore consider
general BMM in the so-called output-sensitive form F ⊙ (AB).
 
  

