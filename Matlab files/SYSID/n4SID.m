%% Setup
clear; clc; close all

%% Load and initialise experiment data
%load experiment data
load('..\..\Data\blockResponse5.mat');

% init output
y1 = t1.long;
y2 = t2.long;
y = [y1; y2];

% init input
u = [u1; u2];

% init time
time = time.long;

%% create Hankel matrices

% Hankel matrix parameters
Ly = size(y,1);
Lu = size(u,1);
s = 150;
n = 6;
N = size(y,2) - 2*s + 1;

% Hankel matrices
Y0sN = zeros(Ly*s,N);
YssN = zeros(Ly*s,N);
U0sN = zeros(Lu*s,N);
UssN = zeros(Lu*s,N);

for i = 1:s
    for j = 1:N 
        Y0sN(Ly*(i-1)+1:Ly*i,j) = y(:,j+i-1);
        YssN(Ly*(i-1)+1:Ly*i,j) = y(:,j+(i+s)-1);
        U0sN(Lu*(i-1)+1:Lu*i,j) = u(:,j+i-1);
        UssN(Lu*(i-1)+1:Lu*i,j) = u(:,j+(i+s)-1);
    end
end

%% N4SID

%RQ factorization;
Zn = [U0sN; Y0sN];
A = [UssN; Zn; YssN];
R = triu(qr(A'))';
% R11 = R(1:Lu*s,1:Lu*s);
R22 = R(Lu*s+1:(2*Lu+Ly)*s, Lu*s+1:Lu*s+(Lu+Ly)*s);
R32 = R((2*Lu+Ly)*s+1:end, Lu*s+1:(2*Lu+Ly)*s);

[~,S,V] = svd(R32/(R22)*Zn,'econ');
S = S(1:n,1:n);
V = V(:,1:n);

XsN = S^(0.5)*V';

%% least-squares - system matrices
% YssNm1 = y(:,s+1:N-1+s);
Ys1Nm1 = YssN(1:Ly,1:end-1);
% UssNm1 = u(:,s+1:N-1+s);
Us1Nm1 = UssN(1:Lu,1:end-1);


x0 = [XsN(:,1:end-1); Us1Nm1];
x1 = [XsN(:,2:end); Ys1Nm1];

F = x0'\x0'*x1/x0;
At = F(1:n,1:n);
Bt = F(1:n,n+1:end);
Ct = F(n+1:end,1:n);
Dt = F(n+1:end,n+1:end);

%% Create the system
xs = zeros(n,numel(time));

xs(:,1) = XsN(:,1);
for i = 1:numel(time)-1
   xs(:,i+1) = At*xs(:,i) + Bt*u(:,i);
end

ys = Ct*xs+Dt*u;


%%
% Compute residuals W and V
W = XsN(:,2:end)-At*XsN(:,1:end-1) - Bt*Us1Nm1;
V = Ys1Nm1-Ct*XsN(:,1:end-1)-Dt*Us1Nm1;

% Covariance of W and V
Q = 1/N*(W*W');
R = 1/N*(V*V');
S = 1/N*(W*V');
% Find the kalman gain by solving the Riccati equation
[~,~,Ks] = dare(At',Ct',Q,R,S);
Ks=Ks';

%% Save system matrices in a struct
NSID.At = At;
NSID.Bt = Bt;
NSID.Ct = Ct;
NSID.Dt = Dt;
NSID.Ks = Ks;
NSID.ys = ys;
NSID.us = u;

save('..\..\Data\NSID.mat','NSID');

%% Plotting

compN4SID = figure('Name','N4SID Comparison');

subplot(2,1,1)
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$',  'Interpreter', 'Latex')
title('Temperature Heater 1')
hold on
plot(time,y1,'b.')
plot(time,ys(1,:),'r','Linewidth',2)
hold off
legend({'y_{1,exp}','y_{1,N4SID}'},'Location','northeast')

subplot(2,1,2)
xlabel('Time in s')
ylabel('Temperature in $^{\circ}C$', 'Interpreter', 'Latex')
title('Temperature Heater 2')
hold on
plot(time,y2,'b.')
plot(time,ys(2,:),'r','Linewidth',2)
hold off
legend({'y_{2,exp}','y_{2,N4SID}'},'Location','northeast')

sgtitle('Comperison experiment data and N4SID output');

% save figure
compN4SID = gcf;
compN4SID.Renderer = 'painters';
saveas(compN4SID, '..\..\Latex\images\SYSID\compN4SID', 'svg');

% Singular values
singularValN4SID = figure('Name','Singular Values R_{32}R_{22}^-1ZN');
hold on
xlabel('singular values')
ylabel('Amplitude')
title('Singular Values of $R_{32}R_{22}^{-1} Z_N$', 'Interpreter', 'Latex')
plot(svd(R32/(R22)*Zn,'econ'),'ro','LineWidth',1);
set(gca,'YScale','log')
hold off

% save figure
singularValN4SID = gcf;
singularValN4SID.Renderer = 'painters';
saveas(singularValN4SID, '..\..\Latex\images\SYSID\singularValN4SID', 'svg');
