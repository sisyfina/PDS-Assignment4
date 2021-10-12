# Exercise 4 - Parallel & Distributed Computer Systems

Sparse Boolean Matrix Multiplication with two levels of parallelization

In this assignment we will implement Boolean matrix multiplication (BMM) of sparse CSC matrices using at least two forms of parallelism that combined produce the correct result faster than either form of parallelism independently. You may combine shared and distributed memory parallelism either OpenMP or Cilk, together with MPI, and or GPU programming.

We denote the logical conjunction (and), disjunction (or) and negation (not) by , ∨ and , respectively. Given matrix with row index set and column index set , matrix with index sets and . Let , and . Let be the Boolean product
:
Element > c(i,j) = 1 if and only if there exists in the search domain such that , i.e., there is a
match between the binary vector in row of and the binary vector in column of .
If provided a filter on the product , then only those elements with are calculated. The filtered
product can be expressed as , where denotes the Hadamard product, i.e., element-wise product. The
absence of a prescribed filter is equivalent to the case where has no zero elements. We may therefore consider
general BMM in the so-called output-sensitive form .


∧ ∨ ¬
A = A(I, K) = [a(i, k)] I K B = B(K, J) = [b(k, j)]
K J ni = |I| nj = |J| nk = |K| C = C(I, J) = [c(i, j)]
C = F ⊙ (AB) c(i, j) = ⋁k∈K a(i, k) ∧ b(k, j), k ∈ K, i ∈ I, j ∈ J, F(i, j) ≠ 0
