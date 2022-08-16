function [theta_hat, w_hat] = tensor_estimator(alpha, beta, gamma, N)
% TENSOR_ESTIMATOR  Apply a spectral method of moments
% estimator for an image with N regions with label distribuitions
%           
%   PARAMS:
%   - alpha, beta, gamma: image first, second and third moments
%   - N: Number Of regions
%   - l: number of labels in the image
%
%   RETURNS:
%   - theta_hat: N-by-l matrix with the prob distribuition estimator for
%   each region    
%   - w_hat: l dimensional with the normalized image region sizes.

k = N;

n = size(gamma, 3);
[U, S, ~] = svd(beta);
E = U(:, 1:k) * sqrt(S(1:k, 1:k));
D = pinv(E);

H = zeros(k, k, n);
for r = 1:n
    H(:, :, r) = D * gamma(:, :, r) * D';
end
O = compute_O(H, N);


theta_hat = E * O;
theta_hat = bsxfun(@rdivide, abs(theta_hat), sum(abs(theta_hat)))';

w_hat = pinv(theta_hat') * alpha;
w_hat = bsxfun(@rdivide, abs(w_hat), sum(abs(w_hat)));
end

function O = compute_O(H, N)
    if N == 2
        h = H(1, 2, :); f = H(1, 1, :) - H(2, 2, :);

        c1 = sum(4 * h.^2 - f.^2);  c2 = -4 * sum(f .* h);
        c3 = 2 * sum(f .* h);       c4 = sum(f .^ 2) - 4 * sum(h .^ 2);
        c5 = sum(h .^ 2);

        x = -1:1e-5:1;
        y  = x.^4 * c1 + x.^3.* sqrt(1 - x.^ 2) * c2 + ...
             x .* sqrt(1 - x.^2) * c3 + x.^2 * c4 + c5;
        a = x(y == min(y));

        O = [sqrt(1 - a^2), a; -a, sqrt(1 - a^2)];
    else
        [O, ~] = rjd(reshape(H, N, N * size(H, 3)));
    end
end

function [V, qDs] = rjd(A, threshold)
%***************************************
% joint diagonalization (possibly
% approximate) of REAL matrices.
%***************************************
% This function minimizes a joint diagonality criterion
% through n matrices of size m by m.
%
% Input :
% * the  n by nm matrix A is the concatenation of m matrices
%   with size n by n. We denote A = [ A1 A2 .... An ]
% * threshold is an optional small number (typically = 1.0e-8 see below).
%
% Output :
% * V is a n by n orthogonal matrix.
% * qDs is the concatenation of (quasi)diagonal n by n matrices:
%   qDs = [ D1 D2 ... Dn ] where A1 = V*D1*V' ,..., An =V*Dn*V'.
%
% The algorithm finds an orthogonal matrix V
% such that the matrices D1,...,Dn  are as diagonal as possible,
% providing a kind of `average eigen-structure' shared
% by the matrices A1 ,..., An.
% If the matrices A1,...,An do have an exact common eigen-structure
% ie a common othonormal set eigenvectors, then the algorithm finds it.
% The eigenvectors THEN are the column vectors of V
% and D1, ...,Dn are diagonal matrices.
%
% The algorithm implements a properly extended Jacobi algorithm.
% The algorithm stops when all the Givens rotations in a sweep
% have sines smaller than 'threshold'.
% In many applications, the notion of approximate joint diagonalization
% is ad hoc and very small values of threshold do not make sense
% because the diagonality criterion itself is ad hoc.
% Hence, it is often not necessary to push the accuracy of
% the rotation matrix V to the machine precision.
% It is defaulted here to the square root of the machine precision.
%
%
% Author : Jean-Francois Cardoso. cardoso@sig.enst.fr
% This software is for non commercial use only.
% It is freeware but not in the public domain.
% A version for the complex case is available
% upon request at cardoso@sig.enst.fr
%-----------------------------------------------------
% Two References:
%
% The algorithm is explained in:
%
%@article{SC-siam,
%   HTML =	"ftp://sig.enst.fr/pub/jfc/Papers/siam_note.ps.gz",
%   author = "Jean-Fran\c{c}ois Cardoso and Antoine Souloumiac",
%   journal = "{SIAM} J. Mat. Anal. Appl.",
%   title = "Jacobi angles for simultaneous diagonalization",
%   pages = "161--164",
%   volume = "17",
%   number = "1",
%   month = jan,
%   year = {1995}}
%
%  The perturbation analysis is described in
%
%@techreport{PertDJ,
%   author = "{J.F. Cardoso}",
%   HTML =	"ftp://sig.enst.fr/pub/jfc/Papers/joint_diag_pert_an.ps",
%   institution = "T\'{e}l\'{e}com {P}aris",
%   title = "Perturbation of joint diagonalizers. Ref\# 94D027",
%   year = "1994" }
%
%
%
[m,nm] = size(A);
V = eye(m);

if nargin==1, threshold = sqrt(eps); end

encore=1;
while encore, encore=0;
    for p = 1:m-1
        for q = p+1:m
            
            % Computation of G (rotations)
            g = [ A(p,p:m:nm)-A(q,q:m:nm) ; A(p,q:m:nm)+A(q,p:m:nm) ];
            g = g * g';
            ton = g(1,1)-g(2,2); toff = g(1,2) + g(2,1);
            theta = 0.5 * atan2( toff , ton+sqrt(ton*ton+toff*toff) );
            c = cos(theta); s = sin(theta);
            encore = encore | (abs(s) > threshold);
            
            % Update of the A and V matrices
            if (abs(s) > threshold)
                Mp = A(:,p:m:nm); Mq = A(:,q:m:nm);
                A(:,p:m:nm) = c * Mp + s * Mq; A(:,q:m:nm) = c * Mq - s * Mp;
                rowp = A(p,:); rowq = A(q,:);
                A(p,:) = c * rowp + s * rowq; A(q,:) = c * rowq - s * rowp;
                temp = V(:,p); V(:,p) = c * V(:,p) + s * V(:,q); V(:,q) = c * V(:,q) - s * temp;
            end
        end
    end
end
qDs = A;

end