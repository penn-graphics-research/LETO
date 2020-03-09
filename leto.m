function leto(nelx, nely, volfrac, penal)
%% Material properties
E0 = 1;
Emin = 1e-6;
nu = 0.3;
%% Initialize carrier particles and load setup.
[xc, yc] = meshgrid(-2:0.5:nelx+2, -2:0.5:nely+2);
carriers = [xc(:)'; yc(:)'; 0.16 * volfrac * ones(1, numel(xc))]; % each column corresponds to a single carrier particle.
F = sparse(2 * (nely + 1) * (nelx + 1), 1, -1, 2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1), 1);
fixeddofs = 1:2*nely;
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% Physical model
dPdF = E0 / (1+nu) / (1-2*nu) * [1-nu, 0, 0, nu; 0, 1/2-nu, 1/2-nu, 0; 0, 1/2-nu, 1/2-nu, 0; nu, 0, 0, 1-nu];
KE = {zeros(8), zeros(8), zeros(8), zeros(8)};
localXq = [0 0 1 1; 0 1 0 1] * 0.5 + 0.25;
for iq = 1:4
    xq = localXq(:, iq);
    dw = [-xq(2), xq(2), 1-xq(2), -(1-xq(2)); 1-xq(1), xq(1),  -xq(1), -(1-xq(1))];
    for in = 1:4 % iterate cached nodes
    for jn = 1:4
        dFdX = zeros(2);
        for id = 1:2 % iterate dim
        for jd = 1:2
            dFdX = dFdX + dPdF(2*id-1: 2*id, 2*jd-1: 2*jd) * dw(id, in) * dw(jd, jn);
        end
        end
        KE{iq}(2*in-1:2*in, 2*jn-1:2*jn) = KE{iq}(2*in-1:2*in, 2*jn-1:2*jn) + 0.25 * dFdX;
    end
    end
end
%% Optimization
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
xold1 = carriers; xold2 = xold1; low = xold1; upp = xold1; % mma intermediate vars
iter = 0; changes = 1e10 * ones(10,1); relCh = 1;
while (relCh > 1e-3 && iter < 10)
    iter = iter + 1;
    dmPow = min(ceil(iter/20), 10); % power of the density map
    baseNode = round(carriers(1:2, :)) + 1;
    [pollutedCellXLocal, pollutedCellYLocal] = meshgrid(-3:2, -3:2);
    pollutedCellLocal = [pollutedCellXLocal(:)'; pollutedCellYLocal(:)']; % each carrier will affect 6 * 6 cells
    pollutedCell = repmat(baseNode, 36, 1) + pollutedCellLocal(:);
    pollutedCell = reshape(pollutedCell, 2, []);
    pollutedCell(1, :) = mod(pollutedCell(1, :)-1, nelx)+1;
    pollutedCell(2, :) = mod(pollutedCell(2, :)-1, nely)+1; % correct negative index;
    pollutedCellIndex = (pollutedCell(1, :) - 1) * nely + pollutedCell(2, :);
    pollutedQuadIndex = 4 * (pollutedCellIndex) + [-3;-2;-1;0];
    quadPos = reshape((repmat(pollutedCell, 4, 1) - 1 + localXq(:)), 72 * 4, []);
    carrierShift = reshape(repmat(carriers(1:2,:),36 * 4, 1) - quadPos, 2, []);
    R = vecnorm(carrierShift)/1.25; % the spline radius is set as 1.25 dx
    W = zeros(size(R)); dW_div_R = zeros(size(R));
    seg1 = R < 1; seg2 = (R < 2) & (R >= 1); 
    W(seg1) = 15 / (7 * pi) * (R(seg1).^3 / 2 - R(seg1).^2 + 2 / 3);
    W(seg2) = 5 / (14 * pi) * (2 - R(seg2)).^3;
    dW_div_R(seg1) = 15 / (7 * pi) * (3 * R(seg1) / 2 - 2);
    dW_div_R(seg2) = - 15 / (14 * pi) * (2 - R(seg2)) .^ 2 ./ R(seg2);
    accuRho = reshape(W, 36 * 4, []) .* carriers(3,:);
    quadRhoNaive = accumarray(pollutedQuadIndex(:), accuRho(:));
    [quadRho, dDM] = density_map(quadRhoNaive, dmPow);
    quadRho = reshape(quadRho, 4, nelx * nely); dDM = reshape(dDM, 4, nelx * nely);
    dQuadRhoNaiveCom = [(reshape(repmat(carriers(3, :), 36*4, 1), size(W)) .* dW_div_R / 1.25^2) .* carrierShift; W];   
    sK = zeros(numel(iK),1);
    for iq = 1:4
        sK = sK + reshape(KE{iq}(:)*(Emin+ quadRho(iq,:).^penal * (E0-Emin)),64*nelx*nely,1);
    end
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs); 
    ce = zeros(size(quadRho));
    for iq = 1:4
        ce(iq, :) = sum((U(edofMat)*KE{iq}).*U(edofMat), 2);
    end
    c = sum(ce .* (Emin+quadRho.^penal*(E0-Emin)) , [1,2]);
    dCdRhoNaive = -penal*(E0-Emin) * reshape(quadRho.^(penal-1).*dDM.*ce,[], 1);
    dC = reshape(sum(reshape(dQuadRhoNaiveCom .* dCdRhoNaive(pollutedQuadIndex(:))', 3, 36*4, []), 2), [], 1);
    v = sum(quadRho, [1,2])/numel(quadRho);
    dVdRhoNaive = reshape(dDM ./ numel(quadRho), [], 1);
    dV = reshape(sum(reshape(dQuadRhoNaiveCom .* dVdRhoNaive(pollutedQuadIndex(:))', 3, 36*4, []), 2), [], 1);
    stepPos = min(nelx, nely) / 40;
    xmin = [carriers(1:2, :) - stepPos; max(carriers(3,:) - 0.5, 0)];
    xmax = [carriers(1:2, :) + stepPos; carriers(3,:) + 0.5];
    dC_norm = max(abs(dC(:))); dV_norm = max(abs(dV(:)));
    [carriers(:), low(:), upp(:), xold1(:), xold2(:)] = mmaUpdate(carriers(:), dC(:)/dC_norm, ... 
        (v - volfrac)/dV_norm, dV(:)/dV_norm, xmin(:), xmax(:), iter, low(:), upp(:), xold1(:), xold2(:));
    changes(mod(iter-1, 10) + 1) = c; relCh = min(1,(max(changes) - min(changes)) / min(changes));
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f dC_norm: %7.3f relCh: %7.3f\n',iter,c,v, dC_norm, relCh);
    %% Visualization
    rho_visual = zeros(2 * nely, 2 * nelx);
    rho_visual(1:2:2*nely, 1:2:2*nelx) = reshape(quadRho(1, :), nely, nelx);
    rho_visual(2:2:2*nely, 1:2:2*nelx) = reshape(quadRho(2, :), nely, nelx);
    rho_visual(1:2:2*nely, 2:2:2*nelx) = reshape(quadRho(3, :), nely, nelx);
    rho_visual(2:2:2*nely, 2:2:2*nelx) = reshape(quadRho(4, :), nely, nelx);
    colormap(gray); imagesc(1-rho_visual); caxis([0 1]); axis equal; axis off; drawnow;
end
end
function [rho, drho] = density_map(x, k)
    rho = ones(length(x), 1); drho = zeros(length(x), 1);
    if k == 1
        seg1 = x <= 0.9; seg2 = (0.9 < x) & (x < 1.1);
        rho(seg1) = x(seg1); rho(seg2) = -2.5 * (x(seg2) - 0.9).^ 2 + x(seg2);
        drho(seg1) = 1; drho(seg2) = -5 * (x(seg2) - 0.9) + 1; 
    else
        seg1 = x <= 0.5; seg2 = (x > 0.5) & (x < 1);
        rho(seg1) = 0.5 * (2 * x(seg1)) .^ k; rho(seg2) = 1 - 0.5 * (2 - 2 * x(seg2)) .^ k;
        drho(seg1) = k * (2 * x(seg1)) .^ (k-1); drho(seg2) = k * (2 - 2 * x(seg2)) .^ (k-1);
    end
end