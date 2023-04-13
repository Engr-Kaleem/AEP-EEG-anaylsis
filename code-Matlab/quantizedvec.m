function y = quantizedvec(x, q)
x = x(:);
nx = length(x);
nq = length(q);
y = sum(repmat(x,1,nq)>repmat(q,nx,1),2);