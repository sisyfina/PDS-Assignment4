%{
  Compressed Sparse Column (CSC) Format

  #integer                   :: m      ! number of rows (unsymmetric only)
  #integer                   :: n      ! number of columns
  #integer, size(n+1)        :: col    ! column pointers (may have type long)
  #integer, size(ptr(n+1)-1) :: row    ! row indices
  #real,    size(ptr(n+1)-1) :: val    ! numerical values
%}  

function [m, n, row_ind, col_start] = cscvarC(a)

[m, n] = size(a);
%em = ones(1,m);
[row,col] = find(a); 
 tmp = size(row);
 noe = tmp(1,1);

 col_start1 = sum(a, 1)';
 row_ind = row - ones(noe,1);
 col_start = zeros(n+1,1);

 for i=1:1:n
     col_start(i+1,1) = col_start1(i,1) + col_start(i,1);
 end


% open your file for writing
 fid = fopen('row_indC.txt','wt');
 % write the matrix
 %myData = randi(255,10,3);
 if fid > 0
     fprintf(fid,'%d,',row_ind');
     fclose(fid);
 end
%
 fid = fopen('col_startC.txt','wt');
 if fid > 0
     fprintf(fid,'%d,',col_start');
     fclose(fid);
 end
     
 fid = fopen('mC.txt','wt');
 if fid > 0
     fprintf(fid,'%d',m');
     fclose(fid);
 end
 
 fid = fopen('nC.txt','wt');
 if fid > 0
     fprintf(fid,'%d',n');
     fclose(fid);
 end


end
    