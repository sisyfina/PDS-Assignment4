%{
% MATLAB uses ordinary multiplication with doubles
% To emulate BMM, we compare against 0 to get true/false values
% do not implement your code with doubles!
function C = bmm(A,B)
C = (A*B)>0;
end

function C = bmmfiltered(A,B,F)
C = ( F.*(A*B) )>0;
end

%}

% parallila

n = 2e6;
%n = 3;
d = 2; % approximate number of true elements per row

A = sprand( n, n, d/n ) > 0;
cscvarA(A);
B = sprand( n, n, d/n ) > 0;
cscvarB(B);
%write data to file????????????????

tic;
C = (double(A)*double(B)) > 0;
cscvarC(C);
toc;