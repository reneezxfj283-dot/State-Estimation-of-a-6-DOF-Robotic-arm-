function X = Gauss(M,S,N)

  if nargin < 3
    N = 1;
  end
  
  L = chol(S)';
  X = repmat(M,1,N) + L*randn(size(M,1),N);