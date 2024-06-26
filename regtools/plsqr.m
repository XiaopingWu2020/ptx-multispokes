function [X,rho,eta,F] = plsqr(A,L,W,b,k,reorth,sm)
%PLSQR "Preconditioned" version of the LSQR Lanczos bidiagonalization algorithm.
%
% [X,rho,eta,F] = plsqr(A,L,W,b,k,reorth,sm)
%
% Performs k steps of the `preconditioned' LSQR Lanczos
% bidiagonalization algorithm applied to the system
%    min || (A*L_p) x - b || ,
% where L_p is the A-weighted generalized inverse of L.  Notice
% that the matrix W holding a basis for the null space of L must
% also be specified.
%
% The routine returns all k solutions, stored as columns of
% the matrix X.  The solution seminorm and the residual norm are
% returned in eta and rho, respectively.
%
% If the generalized singular values sm of (A,L) are also provided,
% then glsqr computes the filter factors associated with each step
% and stores them columnwise in the matrix F.
%
% Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS,
%    reorth = 2 : Householder-reorthogonalization.

% References: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
% sparse linear equations and sparse least squares", ACM Trans.
% Math. Software 8 (1982), 43-71.
% P. C. Hansen, "Rank-Deficient and Discrete Ill-Posed Problems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM, 05/26/93.

% The fudge threshold is used to prevent filter factors from exploding.
fudge_thr = 1e-4;

% Initialization
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==5), reorth = 0; end
if (nargout==4 & nargin<7), error('Too few input arguments'), end
[m,n] = size(A); X = zeros(n,k); [pp,n1] = size(L);
if (n1 ~= n | m < n | n < pp)
  error('Incorrect dimensions of A and L')
end
if (reorth==0)
  UV = 0;
elseif (reorth==1)
  U = zeros(m,k); V = zeros(pp,k); UV = 1;
elseif (reorth==2)
  if (k>=n), error('No. of iterations must satisfy k < n'), end
  UV = 0; HHU = zeros(m,k); HHV = zeros(pp,k);
  HHalpha = zeros(1,k); HHbeta = HHalpha;
else
  error('Illegal reorth')
end
if (nargout > 1)
  eta = zeros(k,1); rho = eta;
  c2 = -1; s2 = 0; xnorm = 0; z = 0;
end
if (nargin==7)
  [ls,ms] = size(sm);
  F = zeros(ls,k); Fv = zeros(ls,1); Fw = Fv;
  s = (sm(:,1)./sm(:,2)).^2;
end

% Prepare for computations with L_p.
[NAA,x_0] = pinit(W,A,b);

% By subtracting the component A*x_0 from b we insure that
% the corrent residual norms are computed.
b = b - A*x_0;

% Prepare for LSQR iteration.
v = zeros(pp,1); x = v; beta = norm(b);
if (beta==0), error('Right-hand side must be nonzero'), end
if (reorth==2)
  [beta,HHbeta(1),HHU(:,1)] = gen_hh(b);
end
u = b/beta; if (UV), U(:,1) = u; end
r = ltsolve(L,A'*u,W,NAA); alpha = norm(r);
if (reorth==2)
  [alpha,HHalpha(1),HHV(:,1)] = gen_hh(r);
end
v = r/alpha; if (UV), V(:,1) = v; end
phi_bar = beta; rho_bar = alpha; w = v;
if (nargin==7), Fv = s/(alpha*beta); Fw = Fv; end

% Perform Lanczos bidiagonalization with/without reorthogonalization.
for i=2:k+1

  alpha_old = alpha; beta_old = beta;

  % Compute (A*L_p)*v - alpha*u.
  p = A*lsolve(L,v,W,NAA) - alpha*u;
  if (reorth==0)
    beta = norm(p); u = p/beta;
  elseif (reorth==1)
    for j=1:i-1, p = p - (U(:,j)'*p)*U(:,j); end
    beta = norm(p); u = p/beta;
  else
    for j=1:i-1
      p(j:m) = app_hh(p(j:m),HHbeta(j),HHU(j:m,j));
    end
    [beta,HHbeta(i),HHU(i:m,i)] = gen_hh(p(i:m));
    u = zeros(m,1); u(i) = 1;
    for j=i:-1:1
      u(j:m) = app_hh(u(j:m),HHbeta(j),HHU(j:m,j));
    end
  end

  % Compute L_p'*A'*u - beta*v.
  r = ltsolve(L,A'*u,W,NAA) - beta*v;
  if (reorth==0)
    alpha = norm(r); v = r/alpha;
  elseif (reorth==1)
    for j=1:i-1, r = r - (V(:,j)'*r)*V(:,j); end
    alpha = norm(r); v = r/alpha;
  else
    for j=1:i-1
      r(j:pp) = app_hh(r(j:pp),HHalpha(j),HHV(j:pp,j));
    end
    [alpha,HHalpha(i),HHV(i:pp,i)] = gen_hh(r(i:pp));
    v = zeros(pp,1); v(i) = 1;
    for j=i:-1:1
      v(j:pp) = app_hh(v(j:pp),HHalpha(j),HHV(j:pp,j));
    end
  end

  % Store U and V if necessary.
  if (UV), U(:,i) = u; V(:,i) = v; end

  % Construct and apply orthogonal transformation.
  rrho = pythag(rho_bar,beta); c1 = rho_bar/rrho;
  s1 = beta/rrho; theta = s1*alpha; rho_bar = -c1*alpha;
  phi = c1*phi_bar; phi_bar = s1*phi_bar;

  % Compute solution norm and residual norm if necessary;
  if (nargout > 1)
    delta = s2*rrho; gamma_bar = -c2*rrho; rhs = phi - delta*z;
    z_bar = rhs/gamma_bar; eta(i-1) = pythag(xnorm,z_bar);
    gamma = pythag(gamma_bar,theta);
    c2 = gamma_bar/gamma; s2 = theta/gamma;
    z = rhs/gamma; xnorm = pythag(xnorm,z);
    rho(i-1) = abs(phi_bar);
  end

  % If required, compute the filter factors.
  if (nargin==7)

    if (i==2)
      Fv_old = Fv;
      Fv  = Fv.*(s - beta^2 - alpha_old^2)/(alpha*beta);
      F(:,i-1) = (phi/rrho)*Fw;
    else
      tmp = Fv;
      Fv = (Fv.*(s - beta^2 - alpha_old^2) - ...
                 Fv_old*alpha_old*beta_old)/(alpha*beta);
      Fv_old = tmp;
      F(:,i-1) = F(:,i-2) + (phi/rrho)*Fw;
    end
    if (i > 3)
      f = find(abs(F(:,i-2)-1) < fudge_thr & abs(F(:,i-3)-1) < fudge_thr);
      if (length(f) > 0), F(f,i-1) = ones(length(f),1); end
    end
    Fw = Fv - (theta/rrho)*Fw;

  end

  % Update the solution.
  x = x + (phi/rrho)*w; w = v - (theta/rrho)*w;
  X(:,i-1) = lsolve(L,x,W,NAA) + x_0;

end
