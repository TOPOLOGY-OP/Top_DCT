function Top_DCT(nelx,nely,volfra)
%% SETUP PARAMETERS
E0 = 200.0e5; nu = 0.3; 
tolne = nelx*nely; tolnd = (nelx+1)*(nely+1); 
tolvol = tolne;
%% PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:tolnd,1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,tolne,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1], tolne,1);
iK = reshape(kron(edofMat,ones(8,1))', 64*tolne,1);
jK = reshape(kron(edofMat,ones(1,8))', 64*tolne,1);
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE0 = 1/(1-nu^2)/24*([A11 A12; A12' A11]+nu*[B11 B12; B12' B11]);
% SETUP BOUNDARY CONDITIONS
F = sparse(2*(nelx*(nely+1)+1), 1, -1000, 2*tolnd, 1);
U = zeros(2*tolnd,1);
fixeddofs = [1*2*(nely+1),2*(nelx)*(nely+1)+1:2:2*(nelx+1)*(nely+1)];
freedofs = setdiff(1:2*tolnd,fixeddofs);
%% INITIALIZE ITERATION
x = repmat(volfra,nely,nelx);
loop = 0; obj = 0.;penal = 3;beta = 0; Cs = [];  volf = [];
%% Call the DCT decomposition
locations = []; c = 40;
for i = 1 : c
    for j = 1 : (i/(nelx/nely))
        temp = (c + 1 - i - 1) * nely + j;
        locations = cat(1,locations,temp);
    end
end
Q = dct(x,[],1);
R = dct(Q,[],2);
X = R(:);
[~,ind] = sort(abs(R(:)),'descend');
coeffs = 1;
while norm(X(ind(1:coeffs)))/norm(X) < 0.999
   coeffs = coeffs + 1;
end
coeffs = length(locations);
R(abs(R) < X(coeffs)) = 0;
disp(["design nums is " + num2str(length(find(R~=0)))]);
change = 1.; ichange = 1; n = coeffs;
x_designs = R(locations(1:coeffs));xmin = -1000*ones(n,1); xmax = 1000*ones(n,1);
low = xmin; upp = xmax;
xold1 = x_designs;  xold2 = x_designs; clf;
%% DIFFERENCE
S = idct(R,[],2); T = idct(S,[],1);
eIntopMat = zeros(coeffs,tolne);delta = 1e-6;
for i = 1: n
    R_temp = R(:);
    R_temp(locations(i)) = R_temp(locations(i)) + delta;
    S_temp = idct(reshape(R_temp, nely, nelx),[],2);
    T_temp = idct(S_temp,[],1);
    eIntopMat(i,:) = (T_temp(:) - T(:))/delta;
end
%% START ITERATION
while (change>=0.005 || beta<40) && loop < 250
    loop = loop + 1;
    objold = obj;
    ePhi = eIntopMat'*x_designs(:);
    [ePhiProj, etha] = THRESHOLD(ePhi, beta);
    edproj = DERIVATIVE_OF_THRESHOLD(ePhi, beta, etha);
    %% DO FINITE ELEMENT ANALYSIS
    sK = reshape(KE0(:)*(max(1e-9,ePhiProj(:)'.^penal)*E0), 64*tolne, 1);
    K = sparse(iK,jK,sK); K=(K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %% EVALUATE OPTIMIZATION RESPONSES AND DO SENSITIVITY ANALYSIS
    obj = F'*U;
    ce = sum((U(edofMat)*KE0).*U(edofMat), 2);
    dcdx = eIntopMat*(-0.5*penal*E0*(ePhiProj).^(penal-1).*ce.*edproj);
    vol = sum(ePhiProj(:));
    voldgdx = eIntopMat*edproj;
    %% UPDATE DESIGN VARIABLES
    m = 1;
    cc = 10000*ones(m,1); d = zeros(m,1); a0 = 1; a = zeros(m,1);
    fval = zeros(m, 1);
    fval(1) = 100*(vol/tolvol-volfra);
    dfdx = zeros(m, n);
    dfdx(1,:) = 100*voldgdx/tolvol;
    [xmma,~,~,~,~,~,~,~,~,low,upp]=...
        mmasub(m,n,loop,x_designs(:),xmin,xmax,xold1,xold2, ...
        obj,dcdx,fval,dfdx,low,upp,a0,a,cc,d);
    xold2 = xold1(:); xold1 = x_designs(:); x_designs = xmma(:);
    %% TUNE PROJECTION PARAMETER
    change = abs(obj-objold)/obj;
    if change < 0.005 && loop > 30
        ichange = ichange+1;
    else
        ichange = 1;
    end
    if mod(ichange,3)==0
        beta = min(beta+1.5,40);
    end
    %% PRINT RESULTS AND PLOT DENSITIES
    fprintf([' It.:%5i Obj.:%9.4f Vol:%7.4f numdesvars :%5i' ...
        ' beta:%5.1f ch.:%6.3f\n'],...
        loop,obj,vol/tolvol,coeffs,beta,change);
    figure(1); 
    displayx = zeros(nely, 2*nelx);
    displayx(:, 1:nelx) = reshape(ePhiProj, nely, nelx);
    displayx(:, nelx+1:end) = displayx(:, nelx:-1:1);
    colormap(gray); clims=[-1 0];imagesc(-displayx,clims); axis equal; axis tight;
    title('Elemental density distribution');
    set(gca,'XTick',[0 1e5]);set(gca,'YTick',[0 1e5]); 
    Cs = cat(2, Cs, obj);
    volf = cat(2, volf, vol/tolvol);
    plotConvergence(1:loop, Cs, volf, 'c');
    pause(1e-6);
end
end
