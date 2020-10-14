%% 1. Parameters
gn = (1.015)^(1/4)-1;
gz = (1.016)^(1/4)-1;
beta = .9722^(1/4);
delta = 1-(1-.0464)^(1/4);
psi = 2.24;
sigma = 1.000001;
theta = .35;

P = [0.98, -0.0138, -0.0117, 0.192; ...
-0.033, 0.956, -0.0451, 0.0569; ...
-0.0702, -0.0460, 0.896, 0.104;...
0.00481 -0.00811, 0.0488, 0.971];

sc = 1.05; % scale of the P matrix to make the AR(1) stationary
P = P./sc;
Sbar = [-0.0239; 0.328; 0.483;-1.53];
P0 = (eye(4)-P)*Sbar;
Q = [0.0116,0,0,0;0.00141,0.00644,0,0;-0.0105,0.00103,0.0158,0; ...
-0.000575,0.00611,0.0142,0.00458];

param = [gn;gz;beta;delta;psi;sigma;theta];

%% 2. Steady State
zs = exp(Sbar(1));
tauls = Sbar(2);
tauxs = Sbar(3);
gs = exp(Sbar(4));
beth = beta*(1+gz)^(-sigma);

kls = ((1+tauxs)*(1-beth*(1-delta))/(beth*theta))^(1/(theta-1))*zs;
A = (zs/kls)^(1-theta)-(1+gz)*(1+gn)+1-delta; %original
B = (1-tauls)*(1-theta)*kls^theta*zs^(1-theta)/psi; %original
ks = (B+gs)/(A+B/kls); %original
cs = A*ks-gs; %original
ls = ks/kls; %original
ys = ks^theta*(zs*ls)^(1-theta);
xs = ys-cs-gs;

display(['Steady state k/y is ',num2str(ks/ys)]);
display(['Steady state c/y is ',num2str(cs/ys)]);
display(['Steady state hours is ',num2str(ls)]);

%% 3. Call subroutine with residuals:
Z = [log(ks);log(ks);log(ks);log(zs);log(zs);tauls;tauls;
tauxs;tauxs;log(gs);log(gs)];
del = max(abs(Z)*1e-5,1e-8);
for i=1:11;
Zp = Z;
Zm = Z;
Zp(i) = Z(i)+del(i);
Zm(i) = Z(i)-del(i);
dR(i,1) = (res_wedge(Zp,param)-res_wedge(Zm,param))/(2*del(i));
end

%% 4. Solution: log k[t+1] = gamma0 + gammak* log k[t] + gamma* S[t]

a0 = dR(1);
a1 = dR(2);
a2 = dR(3);
b0 = dR(4:2:11)';
b1 = dR(5:2:11)';
tem = roots([a0,a1,a2]);
gammak = tem(find(abs(tem)<1));
gamma = -((a0*gammak+a1)*eye(4)+a0*P')\(b0*P+b1)';
gamma0 = (1-gammak)*log(ks)-gamma'*[log(zs);tauls;tauxs;log(gs)];
Gamma = [gammak;gamma;gamma0];

%% 5. Generate a time series for the 4-dimensional vector of shocks
% according to the specified process, the time series for capital
% according to the above derived policy function, for 100 and 1000
% observations
T = 1000;
t = 100;
st = zeros(T,4);
lkt = zeros(T,1);
lkt(1) = log(ks);
st(1,:) = [Sbar(1) Sbar(2) Sbar(3) Sbar(4)];
sq = 100; %scale the Q matrix
et = randn(T,4);
for i=2:T
st(i,:) = P0+P*st(i-1,:)'+sq*Q*Q'*et(i-1,:)';
lkt(i) = gamma0 + gammak*lkt(i-1)+gamma'*st(i,:)';
end

%% 6. Plotting
subplot(2,2,1)
plot((1:t)./4,exp(lkt(1:t))./ks)
title(['log(kt) for n=',num2str(t)])
xlim([0 t(end)/4]);
subplot(2,2,2)
plot((1:t)./4,st(1:t,:))
xlim([0 t(end)/4]);
title(['log(st) for n=',num2str(t)])
subplot(2,2,3)
plot((1:T)./4,lkt(1:T)./log(ks))
title(['log(kt) for n=',num2str(T)])
xlim([0 T(end)/4]);
subplot(2,2,4)
plot((1:T)./4,st)
xlim([0 t(end)/4]);
title(['log(st) for n=',num2str(T)])
