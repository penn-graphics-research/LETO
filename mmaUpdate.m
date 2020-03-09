function [xval, low, upp, xold1, xold2] = mmaUpdate(xval, dfdx, gx, dgdx, xmin, xmax, iter, low, upp, xold1, xold2)
%% parameters
n = length(xval);
m = length(gx);
asyminit = 0.02;
asymdec = 0.65;
asyminc = 1.05;
xmamieps = 1.0e-5;
epsimin = sqrt(n + m) * 1e-9;
raa0 = 0.0001;
albefa = 0.4;
a = zeros(1,m);
c = 1000 * ones(1,m);

%% Generate the subproblem
if iter < 3
    low = xval - asyminit * (xmax - xmin);
    upp = xval + asyminit * (xmax - xmin);
else
    zzz = (xval - xold1) .* (xold1 - xold2);
    gamma = ones(n, 1);
    gamma(zzz < 0) = asymdec;
    gamma(zzz > 0) = asyminc;
    low = xval - gamma .* (xold1 - low);
    upp = xval + gamma .* (upp - xold1);
    xmami = max(xmamieps, xmax - xmin);
    low = max(low, xval - 10.0 * xmami);
    low = min(low, xval - 0.01 * xmami);
    upp = max(upp, xval + 0.01 * xmami);
    upp = min(upp, xval + 10.0 * xmami);
end

% Set bounds and the coefficients for the approximation
alpha = max(xmin, low + albefa * (xval - low));
beta = min(xmax, upp - albefa * (upp - xval));

% objective function
dfdxp = max(0.0, dfdx); %p0
dfdxm = max(0.0, -1.0 * dfdx); %q0
xmamiinv = 1.0 ./ max(xmamieps, xmax - xmin);
pq = 0.001 * (dfdxp + dfdxm) + raa0 * xmamiinv;
p0 = (upp - xval) .^ 2.0 .* (dfdxp + pq);
q0 = (xval - low) .^ 2.0 .* (dfdxm + pq);

% Constraints
dgdxp = max(0.0, dgdx);
dgdxm = max(0.0, -1.0 * dgdx);
xmamiinv = 1.0 ./ max(xmamieps, xmax - xmin);
pq = 0.001 * (dgdxp + dgdxm) + raa0 * xmamiinv;
pij = (upp - xval).^2 .* (dgdxp + pq);
qij = (xval - low).^2 .* (dgdxm + pq);

b = -gx + sum(pij ./ (upp - xval) + qij ./ (xval - low), 1);
xold2 = xold1;
xold1 = xval;

%% Solve the dual with an interior point method
lam = c/2;
mu = ones(1,m);
tol = epsimin;
epsi = 1.0;
err = 1.0;
while epsi > tol
    loop = 0;
    while (err > 0.9 * epsi) && (loop < 300)
        loop = loop + 1;
        % Set up Newton system 
        % XYZofLAMBDA
        lam = max(0, lam);
        y = max(0, lam - c); % Note y=(lam-c)/d - however d is fixed at one !!
        lamai = dot(lam, a);
        z = max(0, 10 * (lamai - 1)); % SINCE a0 = 1.0
        pjlam = p0 + sum(pij .* lam, 2);
        qjlam = q0 + sum(qij .* lam, 2);
        xval = (sqrt(pjlam) .* low + sqrt(qjlam) .* upp) ./ (sqrt(pjlam) + sqrt(qjlam));
        xval = max(alpha, xval);
        xval = min(beta, xval);
        % DualGrad
        grad = -b - a * z - y + sum(pij ./ (upp - xval) + qij ./ (xval - low), 1);
        grad = - grad - epsi ./ lam;
        % DualHess
        PQ = pij ./ (upp - xval).^2 - qij ./ (xval - low).^2;
        df2 = -1 ./ (2 * pjlam ./ (upp - xval).^3 + 2 * qjlam ./ (xval - low) .^3);
        xp = (sqrt(pjlam) .* low + sqrt(qjlam) .* upp) ./ (sqrt(pjlam) + sqrt(qjlam));
        df2(xp < alpha) = 0;
        df2(xp > beta) = 0;
        % Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
        hess = PQ' * (df2 .* PQ);
        lam = max(0, lam);
        lamai = lamai + dot(lam, a);
        diag_mod = zeros(1, m);
        diag_mod(lam > c) = -1;
        diag_mod = diag_mod - mu./lam;
        hess = hess + diag(diag_mod);
        if lamai > 0
            hess = hess - 10 *  a' * a;
        end
        HessCorr = 1e-4 * trace(hess) / m;
        if HessCorr > -1.0e-7
            HessCorr = -1.0e-7;
        end
        hess = hess + eye(m) * HessCorr;
        % Solve Newton system
        s = hess \ grad';
        % Get the full search direction
        s = [s'; -mu + epsi ./ lam - s' .* mu ./ lam];
        % DualLineSearch
        theta = 1 / max([1.005, -1.01 * s(1, :) ./ lam, -1.01 * s(2, :) ./ mu]);
        lam = lam + theta * s(1,:);
        mu = mu + theta * s(2, :);
        % XYZofLAMBDA
        lam = max(0, lam);
        y = max(0, lam - c); % Note y=(lam-c)/d - however d is fixed at one !!
        lamai = dot(lam, a);
        z = max(0, 10 * (lamai - 1)); % SINCE a0 = 1.0
        pjlam = p0 + sum(pij .* lam, 2);
        qjlam = q0 + sum(qij .* lam, 2);
        xval = (sqrt(pjlam) .* low + sqrt(qjlam) .* upp) ./ (sqrt(pjlam) + sqrt(qjlam));
        xval = max(alpha, xval);
        xval = min(beta, xval);
        % DualResidual
        res = [-b - a * z - y + mu + sum(pij ./ (upp - xval) + qij ./ (xval - low), 1); mu .* lam - epsi];
        err = max(abs(res(:)));
    end
    epsi = epsi * 0.1;
end