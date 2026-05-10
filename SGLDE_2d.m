clear;
% SGLD_2d sampling

% fix the random variable
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

% settings
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 10);
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultTextFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 15);   
set(groot, 'DefaultTextFontSize', 15);    

% cut-off function psi
function psi = PSI(x,k)
    y= abs(k*x);
    eps_val = 1e-10;
    psi = (y>=2).*y+(y>1 & y<2) .* ( exp(-1./(max(y-1,0)+eps_val))...
        ./(exp(-1./(max(y-1,0)+eps_val))+exp(-1./(max(2-y,0)+eps_val)))).*y;
end

% KL-distance calculation

function kle = KLE(X1,X2)
    % histogram
    nx = 100; ny = 100; 
    xlin = linspace(min([X1(:,1);X2(:,1)]), max([X1(:,1);X2(:,1)]), nx);
    ylin = linspace(min([X1(:,2);X2(:,2)]), max([X1(:,2);X2(:,2)]), ny);
    [XX,YY] = meshgrid(xlin, ylin);
    grid_points = [XX(:), YY(:)];

    % kernel density
    p = ksdensity(X1, grid_points, 'Bandwidth', 0.1); 
    q = ksdensity(X2, grid_points, 'Bandwidth', 0.1); 

     % convert density -> discrete probability with cell area
    dx = xlin(2) - xlin(1);
    dy = ylin(2) - ylin(1);
    dA = dx * dy;

    p = p * dA;
    q = q * dA;

    % normalize (numerical)
    p = p / sum(p);
    q = q / sum(q);

    eps_val = 1e-12;
    p = p + eps_val;
    q = q + eps_val;

    kle = sum(p .* log(p ./ q));
end


%% initial data

% U = 1/4 |x|^4 -  1/2 |x|^2, gradU = |x|^2 x - x; beta = 1/2;
N = 2e2;
INI = randn(N,2);
T = 30;   Tt = 7030;
A = 3 : 7;
h = 2.^(-A);
href = 2^(-12); nref = 2^12;
alpha = 1/2;
gamma = 1e-1;
KL = zeros(3, max(size(A)));  % 1:TaE, 2: MTE, 3: MTE-RBM  
Dref = zeros(N,Tt-T,2);

 %% long-time estimation & tamed SGLD
xref = INI;
for i = 0 : Tt*nref
    P = sum(xref.^2,2);
    b = (P-1).*xref;
    drift = - b./(1+PSI(href^alpha*sqrt(sum(b.^2,2)),gamma)); 
    xref = xref + drift.* href + sqrt(href).*randn(N,2);
    if(i>T*nref && (mod(i,nref)==0)) 
        Dref(:,i/nref-T,1) = xref(:,1);
        Dref(:,i/nref-T,2) = xref(:,2);
    end
end


for j = 1 : max(size(A))
    x1 = INI; x2 = INI; x3 = INI;
    D1 = zeros(N,Tt-T,2); D2 = zeros(N,Tt-T,2); D3 = zeros(N,Tt-T,2);
    n = 2^A(1,j);
    for i = 0: Tt*n
        dw = sqrt(h(1,j))*randn(N,2);
        % tamed Euler(TaE)
        P = sum(x1.^2,2);
        b = (P-1).*x1;
        drift = - b./(1+h(1,j)^alpha*sqrt(sum(b.^2,2))); 
        x1 = x1 + drift.* h(1,j) + dw;

        % modified tamed Euler(MTE)
        P = sum(x2.^2,2);
        b = (P-1).*x2;
        drift = - b./(1+PSI(h(1,j)^alpha*sqrt(sum(b.^2,2)),gamma)); 
        x2 = x2 + drift.* h(1,j) + dw;

        % MTE-RBM
        P = sum(x3.^2,2);
        a  = rand(1);
        if(a<0.5)
            b = (2.*P).*x3;
        else 
            b = -2.*x3;
        end
            drift = - b./(1+PSI(h(1,j)^alpha*sqrt(sum(b.^2,2)),gamma)); 
            x3 =x3 + drift.* h(1,j) + dw;

        % sampling
        if(i>T*n && (mod(i,n)==0)) 
            k = i/n - T;
            D1(:,k,1) = x1(:,1); D1(:,k,2) = x1(:,2);
            D2(:,k,1) = x2(:,1); D2(:,k,2) = x2(:,2);
            D3(:,k,1) = x3(:,1); D3(:,k,2) = x3(:,2);
        end
    end
    Dref = reshape(Dref, [], 2);
    D1 = reshape(D1, [], 2);
    D2 = reshape(D2, [], 2);
    D3 = reshape(D3, [], 2);

    %% KL distance (use ksdensity)
    KL(1,j) = KLE(Dref,D1);
    KL(2,j) = KLE(Dref,D2);
    KL(3,j) = KLE(Dref,D3);
end

B = (1/2).^ (0:max(size(A))-1) ;


figure;
loglog(h,sum(KL(:,1))/3.*B);
hold on;
loglog(h,sum(KL(:,1))/3.*(B.^2));
loglog(h,KL(1,:),'-ko');
loglog(h,KL(2,:),'-ks');
loglog(h,KL(3,:),'-kd');
xticks(flip(h));
xticklabels({'2^{-7}', '2^{-6}', '2^{-5}', '2^{-4}', '2^{-3}'});
xlabel('h');
ylabel('KL distance');
legend('slope 1', 'slope 2','TaE','MTE','MTE-RBM');
axis tight ;