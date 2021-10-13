# Exercise 4 - Parallel & Distributed Computer Systems

# Sparse Boolean Matrix Multiplication with two levels of parallelization

In this assignment we will implement Boolean matrix multiplication (BMM) of sparse CSC matrices using at least two forms of parallelism that combined produce the correct result faster than either form of parallelism independently. You may combine shared and distributed memory parallelism either OpenMP or Cilk, together with MPI, and or GPU programming.

We denote the logical conjunction (and), disjunction (or) and negation (not) by ∧, ∨ and ¬, respectively. Given matrix A = A(I, K) = [a(i, k)] with row index set I and column index set K, matrix B = B(K, J) = [b(k, j)] with index sets K and J. Let ni = |I|, nj = |J| and nk = |K|. Let C = C(I, J) = [c(i, j)] be the Boolean product C = F ⊙ (AB) c(i, j) = ⋁k∈K a(i, k) ∧ b(k, j), k ∈ K, i ∈ I, j ∈ J, F(i, j) ≠ 0 

Element c(i,j) = 1 if and only if there exists k in the search domain K such that a(i,k)b(k,j) = 1, i.e., there is a match between the binary vector in row i of A and the binary vector in column j of B.

If provided a filter F on the product AB, then only those elements with F(i,j) = 1 are calculated. The filtered product can be expressed as F ⊙ (AB), where ⊙ denotes the Hadamard product, i.e., element-wise product. The absence of a prescribed filter is equivalent to the case where F has no zero elements. We may therefore consider
general BMM in the so-called output-sensitive form F ⊙ (AB).

Similar to algebraic matrix multiplication, BMM can be written in block version. Assume here for convenience that nI = nJ = nK = n. Let b be the block size. Suppose nb = n/b be an integer greater than b. Partition each index set into nb subsets of equal size b. The bipartite matrices A, B and C are accordingly partitioned into nb × nb block matrices with blocks of size b × b .

Denote by Cp,q the (p, q) block of matrix C, 1 ≤ p,q ≤ nb. The matrix blocks of A, B and F are similarly denoted.
Cp,q := 0;
Then, for s := 1, ⋯ , nb
          Cp,q := Cp,q ∨ Fp,q ⊙ ( Ap,s Bs,q )

There are nb^2 blocks in C, and the are nb^3 block products Ap,s Bs,q in total.

Your implementation does not need to be recursive. Select a block size that provides best performance for matrices that have n more than a million.

You should provide two different interfaces, in the same way as the following two MATLAB routines:

```
% MATLAB uses ordinary multiplication with doubles
% To emulate BMM, we compare against 0 to get true/false values
% do not implement your code with doubles!
function C = bmm(A,B)
C = (A*B)>0;
end
function C = bmmfiltered(A,B,F)
C = ( F.*(A*B) )>0;
end
```

Before proceeding to any parallel implementation, verify that your sequential code is correct and its wall-clock
execution time is comparable to the following MATLAB/Octave commands

```
n = 5e6;
d = 2; % approximate number of true elements per row
A = sprand( n, n, d/n ) > 0;
B = sprand( n, n, d/n ) > 0;
tic;
C = (A*B) > 0;
toc
```

# Optional Improvements

## Output sensitive product

There are redundant operations in BMM by the naive calculation, especially when matrices A and B are dense. Once a true is computed in the accumulation for element c(i, j), the rest of the accumulation is redundant and can be skipped. This also integrates with the prescribed filter F.

Cp,q := 0; Xp,q := Fp,q
for s := 1, ⋯ ,nb
Cp,q := Cp,q ∨ Fp,q ⊙ ( Ap,s Bs,q )
Xp,q := Xp,q ∧ (¬Cp,q ).
The remaining block BMMs, during the accumulation, are adaptively filtered, i.e., made output sensitive.

## Method of Four Russians

The block BMMs at the base level, Cp,q := Cp,q ∨ ( Xp,q ⊙ ( Ap,s Bs,q )), may be accelerated by the Four Russians (4R) algorithm. 

# What to submit
 - A 3-page report in PDF format (any pages after the 3rd one will not be taken into account). Describe in a paragraph two applications of BMM in science and engineering in adequate detail. Describe your parallel algorithm design main points, data distribution, communication patterns and blocking
factor choice decisions. Report execution times of your implementations with respect to size n on the sets of matrices we will specify. Use innovative and imformative plots to convey the information and explain the performance you get 
 - Upload the source code on GitHub, BitBucket, Dropbox, Google Drive, etc. and add a link in your report. 
 - Cite the sources of everything you used. 
 - You may work in groups of two. Submit report twice with both collaborator names listed.
  

