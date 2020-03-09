function leto(nelx, nely, volfrac, penal)
E0 = 1; Emin = 1e-6; nu = 0.3;
[xc, yc] = meshgrid(-2:0.5:nelx+2, -2:0.5:nely+2);
carriers = [xc(:)'; yc(:)'; 0.16 * volfrac * ones(1, numel(xc))]; % each column corresponds to a single carrier particle.
F = sparse(2 * (nely + 1) * (nelx + 1), 1, -1, 2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1), 1);
freedofs = setdiff(1:2*(nely+1)*(nelx+1),1:2*nely);
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
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
xold1 = carriers; xold2 = xold1; low = xold1; upp = xold1; % mma intermediate vars
iter = 0; changes = 1e10 * ones(10,1); relCh = 1;
while (relCh > 1e-3 && iter < 200)
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
    W = 15/(7*pi)*((R.^3/2-R.^2+2/3).*(R < 1) + 1/6*(2-R).^3.*   ((R>=1)&(R<2)));
    dW_div_R =  15/(7*pi)*((3*R/2-2).*(R < 1) - 1/2*((2-R).^2./R).*((R>=1)&(R<2)));
    x = accumarray(pollutedQuadIndex(:), reshape(reshape(W, 36*4,[]).*carriers(3,:),[],1));
    if dmPow == 1
        quadRho = x .* (x<=0.9) + (-2.5*(x-0.9).^2+x) .* ((0.9<x)&(x<1.1)) + (x>=1.1);
        dDM =     1 .* (x<=0.9) + ( -5 *(x-0.9)   +1) .* ((0.9<x)&(x<1.1)); 
    else
        quadRho = 0.5*(2*x).^dmPow    .* (x<=0.5) + (1-0.5*(2-2*x).^dmPow)    .* ((x>0.5)&(x < 1)) + (x>=1);
        dDM =   dmPow*(2*x).^(dmPow-1).* (x<=0.5) + (dmPow*(2-2*x).^(dmPow-1)).*((x>0.5)&(x < 1));
    end
    quadRho = reshape(quadRho, 4, nelx * nely); dDM = reshape(dDM, 4, nelx * nely);
    dQuadRhoNaiveCom = [(reshape(repmat(carriers(3, :),36*4,1), size(W)).*dW_div_R/1.25^2).*carrierShift; W];   
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
    c = sum(ce.*(Emin+quadRho.^penal*(E0-Emin)),[1,2]); v = sum(quadRho,[1,2])/numel(quadRho);
    dCdRhoNaive = -penal*(E0-Emin) * reshape(quadRho.^(penal-1).*dDM.*ce,[], 1);
    dC = reshape(sum(reshape(dQuadRhoNaiveCom .* dCdRhoNaive(pollutedQuadIndex(:))', 3, 36*4, []), 2), [], 1);
    dVdRhoNaive = reshape(dDM ./ numel(quadRho), [], 1);
    dV = reshape(sum(reshape(dQuadRhoNaiveCom .* dVdRhoNaive(pollutedQuadIndex(:))', 3, 36*4, []), 2), [], 1);
    xmin = [carriers(1:2, :) - min(nelx, nely) / 40; max(carriers(3,:) - 0.5, 0)];
    xmax = [carriers(1:2, :) + min(nelx, nely) / 40; carriers(3,:) + 0.5];
    dC_norm = max(abs(dC(:))); dV_norm = max(abs(dV(:)));
    [carriers(:), low(:), upp(:), xold1(:), xold2(:)] = mmaUpdate(carriers(:), dC(:)/dC_norm, ... 
        (v - volfrac)/dV_norm, dV(:)/dV_norm, xmin(:), xmax(:), iter, low(:), upp(:), xold1(:), xold2(:));
    changes(mod(iter-1, 10) + 1) = c; relCh = min(1,(max(changes) - min(changes)) / min(changes));
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f dC_norm: %7.3f relCh: %7.3f\n',iter,c,v, dC_norm, relCh);
    rho_visual = zeros(2 * nely, 2 * nelx);
    rho_visual(1:2:2*nely, 1:2:2*nelx) = reshape(quadRho(1, :), nely, nelx);
    rho_visual(2:2:2*nely, 1:2:2*nelx) = reshape(quadRho(2, :), nely, nelx);
    rho_visual(1:2:2*nely, 2:2:2*nelx) = reshape(quadRho(3, :), nely, nelx);
    rho_visual(2:2:2*nely, 2:2:2*nelx) = reshape(quadRho(4, :), nely, nelx);
    colormap(gray); imagesc(1-rho_visual); caxis([0 1]); axis equal; axis off; drawnow;
end