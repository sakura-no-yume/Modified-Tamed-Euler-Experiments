% 1d modified tamed Euler
% b = -x^3-1.875x, sigma = 0.5x
clear;

% fix the random variable
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

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
    % psi = ((y>=1).*y +(y>0.5 & y<1) .* ( y-0.5)*2) ; 
    psi = (y>=2).*y+(y>1 & y<2) .* ( exp(-1./(max(y-1,0)+eps_val))...
        ./(exp(-1./(max(y-1,0)+eps_val))+exp(-1./(max(2-y,0)+eps_val)))).*y;
end

% truncate function pi_Delta for this case
function pi_Delta = PI_DELTA(x,Delta)
    % since b = -x^3-1.875x, sigma = 0.5, bounded by 8x^3 for x >= 1
    % then h(1) = 1, one may choose kappa(\Detla) = \Delta^{-1/12}
    if(norm(x) < Delta^(-1/12)) 
        pi_Delta = x;
    else
        pi_Delta = Delta^(-1/12)*sign(x);
    end
end
    
% time level;
A=5:9;      T=1;
h=2.^(-A);  href=2^(-15);
n=T./h;     nref=T/href;

% initial data
INI = 1;    
alpha= 0.5; gamma= 1;

N = 1e5; % orbit number
S= max(size(A));
SE = zeros(4,S); % 1:TaE, 2:TrE, 3:MTE, 4:MTE_RBM
WE = zeros(12,S);% 4*3
PD = zeros(12,S);% partial data
PDref = zeros(3,1);
drift = zeros(1,1);
diffusion = zeros(1,1);

for i = 1 : N
    Xref = INI;
    X = INI * ones(1,S);
    X1 = X; X2 = X; X3 = X; X4 = X;
    Wref=randn(1,nref)*sqrt(href);

    for j = 1 : nref
         % reference solution, use MTE with slight step size
         b = -Xref^3-1.875*Xref;
         drift = b/(1+PSI(href^alpha*abs(b),gamma));
         Xref = Xref + href.* drift + 0.5*Xref*Wref(1,j);
    end

    for k = 1: S
        W = sum(reshape(Wref,1,nref/n(1,k),[]),2);
        W = (squeeze(W))';
        % WW=reshape(Wref', [nref/n(1,k),n(1,k)]);
        % WW1=sum(WW,1);
        % W=reshape(WW1,[n(1,k),1])';
        for j = 1 : n(1,k)
            % TaE(Sabanis)
            b = -X1(1,k)^3-1.875*X1(1,k);
            drift= b/(1+abs(h(1,k)^alpha*b));
            X1(:,k) = X1(:,k) + h(1,k) .* drift+0.5.*X1(:,k).* W(1,j);

            % TrE(Mao)
            x = PI_DELTA(X2(1,k),h(1,k));
            drift = -x^3-1.875*x;
            diffusion= 0.5.*x;
            X2(:,k) = X2(:,k) + h(1,k) .* drift + diffusion.* W(1,j);

            % MTE(no batch)
            b = -X3(1,k)^3-1.875*X3(1,k);
            drift = b/(1+PSI(h(1,k)^alpha*abs(b),gamma));
            X3(:,k) = X3(:,k) + h(1,k) .* drift +0.5.*X3(:,k).*W(1,j);


            % MTE_RBM
            a  = rand(1);
            if(a<0.5)
                b= - 2*X4(1,k)^3;
                drift = b/(1+PSI(h(1,k)^alpha*abs(b),gamma));
                X4(:,k) = X4(:,k) + h(1,k) .* drift + 0.5.*X4(:,k).* W(1,j);
            else 
                b = -3.75*X4(1,k);
                drift = b/(1+PSI(h(1,k)^alpha*abs(b),gamma));
                X4(:,k) = X4(:,k) + h(1,k) .* drift + 0.5.*X4(:,k).*W(1,j);
            end
        end
        
    
       
    end

    for k = 1 :S
        % strong error
        SE(1,k)=(SE(1,k)*(i-1) + norm(Xref-X1(:,k))^2 )/i;
        SE(2,k)=(SE(2,k)*(i-1) + norm(Xref-X2(:,k))^2 )/i;
        SE(3,k)=(SE(3,k)*(i-1) + norm(Xref-X3(:,k))^2 )/i;
        SE(4,k)=(SE(4,k)*(i-1) + norm(Xref-X4(:,k))^2 )/i;
        % weak error, 
        % test function f(x) = x;
        PD(1,k) = (PD(1,k).*(i-1)+ X1(:,k)  )/i;
        PD(2,k) = (PD(2,k).*(i-1)+ X2(:,k)  )/i;
        PD(3,k) = (PD(3,k).*(i-1)+ X3(:,k) )/i;
        PD(4,k) = (PD(4,k).*(i-1)+ X4(:,k) )/i;
        % test function f(x) = cos(x);
        PD(5,k) = (PD(5,k).*(i-1)+ cos(X1(:,k))) /i;  
        PD(6,k) = (PD(6,k).*(i-1)+ cos(X2(:,k))) /i; 
        PD(7,k) = (PD(7,k).*(i-1)+ cos(X3(:,k))) /i;
        PD(8,k) = (PD(8,k).*(i-1)+ cos(X4(:,k))) /i;
        % test function f(x) = cos(exp(x));
        PD(9,k) = (PD(9,k).*(i-1)+ cos(exp(X1(:,k)))) /i;  
        PD(10,k) = (PD(10,k).*(i-1)+ cos(exp(X2(:,k)))) /i; 
        PD(11,k) = (PD(11,k).*(i-1)+ cos(exp(X3(:,k)))) /i;
        PD(12,k) = (PD(12,k).*(i-1)+ cos(exp(X4(:,k)))) /i;
    end
    PDref(1,1) = (PDref(1,1).*(i-1)+  Xref )/i;
    PDref(2,1) = (PDref(2,1).*(i-1)+ cos(Xref))/i;
    PDref(3,1) = (PDref(3,1).*(i-1)+ cos(exp(Xref)))/i;
end

for i = 1:4
    WE(i,:)=abs(PD(i,:)-PDref(1,1));
    WE(i+4,:)=abs(PD(4+i,:)-PDref(2,1));
    WE(i+8,:)=abs(PD(8+i,:)-PDref(3,1));
end

disp( 'the initial data is');
disp(INI);
% loglog plot

B = (1/2).^ (0:S-1) ;

figure;
loglog(h,sum(sqrt(SE(1:4,1)))/4.*B);
hold on;
loglog(h,sum(sqrt(SE(1:4,1)))/4.*sqrt(B));
H11 = loglog(h,sqrt(SE(1,:)),'-ko');
H12 = loglog(h,sqrt(SE(2,:)),'-kx');
H13 = loglog(h,sqrt(SE(3,:)),'-ks');
H14 = loglog(h,sqrt(SE(4,:)),'-kd');
xticks(flip(h));
xticklabels({'2^{-9}', '2^{-8}', '2^{-7}', '2^{-6}', '2^{-5}'});
axis tight;
legend([H11,H12,H13,H14],{'TaE','TrE','MTE','MTE-RBM'});
xlabel('h');
ylabel('Strong errors');

figure;
loglog(h,sum(WE(1:4,1))/4.*B);
hold on;
loglog(h,sum(WE(1:4,1))/4.*sqrt(B));
H21 = loglog(h,WE(1,:),'-ko');
H22 = loglog(h,WE(2,:),'-kx');
H23 = loglog(h,WE(3,:),'-ks');
H24 = loglog(h,WE(4,:),'-kd');
xticks(flip(h));
xticklabels({'2^{-9}', '2^{-8}', '2^{-7}', '2^{-6}', '2^{-5}'});
axis tight;
xlabel('h');
ylabel('Weak errors');
legend([H21,H22,H23,H24],{'TaE','TrE','MTE','MTE-RBM'});

figure;
loglog(h,sum(WE(5:8,1))/4.*B);
hold on;
loglog(h,sum(WE(5:8,1))/4.*sqrt(B));
H31 = loglog(h,WE(5,:),'-ko');
H32 = loglog(h,WE(6,:),'-kx');
H33 = loglog(h,WE(7,:),'-ks');
H34 = loglog(h,WE(8,:),'-kd');
xticks(flip(h));
xticklabels({'2^{-9}', '2^{-8}', '2^{-7}', '2^{-6}', '2^{-5}'});
axis tight;
xlabel('h');
ylabel('Weak errors');
legend([H31,H32,H33,H34],{'TaE','TrE','MTE','MTE-RBM'});

figure;
loglog(h,sum(WE(9:12,1))/4.*B);
hold on;
loglog(h,sum(WE(9:12,1))/4.*sqrt(B));
H41 = loglog(h,WE(9,:),'-ko');
H42 = loglog(h,WE(10,:),'-kx');
H43 = loglog(h,WE(11,:),'-ks');
H44 = loglog(h,WE(12,:),'-kd');
xticks(flip(h));
xticklabels({'2^{-9}', '2^{-8}', '2^{-7}', '2^{-6}', '2^{-5}'});
axis tight;
xlabel('h');
ylabel('Weak errors');
legend([H41,H42,H43,H44],{'TaE','TrE','MTE','MTE-RBM'});
