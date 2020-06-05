% Subspace identification using general input sequences and MOESP

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
s = 60;
n = 7;
N = size(y,2) - 2*s + 1;

% Hankel matrices
Y0sN = zeros(Ly*s,N);
U0sN = zeros(Lu*s,N);
for i = 1:s
    for j = 1:N 
        Y0sN(Ly*(i-1)+1:Ly*i,j) = y(:,j+i-1);
        U0sN(Lu*(i-1)+1:Lu*i,j) = u(:,j+i-1);
    end
end

figure
plot(svd(Y0sN),'ro','LineWidth',1);
set(gca,'YScale','log')

%% MOESP method (p312)

% RQ factorisation
R = triu(qr([U0sN;Y0sN]'))';
R22 = R(s*Lu+1:end, s*Lu+1:s*Lu+s*Ly);

% R22 SVD
[U, S, V] = svd(R22);
Un = U(:,1:n);

% Determine Ct and At
Ct = Un(1:Ly,:);
At = Un(1:(s-1)*Ly,:)\Un(Ly+1:s*Ly,:);

% Determine Bt and Dt
tic

phi_1 = [];
phi_2 = [];
phi_3 = [];
for k = 0:size(y,2)-1
   phi_1 = [phi_1; Ct*At^(k)];
   
   phi_sum = zeros(Ly,2*n);
   for tau=0:k-1
      phi_sum = phi_sum + [u(1,tau+1).*Ct*At^(k-tau-1), u(2,tau+1).*Ct*At^(k-tau-1)] ;  
   end
   phi_2 = [phi_2; phi_sum];
   
   phi_3 = [phi_3; [u(1,k+1).*eye(Ly), u(1,k+1).*eye(Ly)]];
end

phi = [phi_1 phi_2 phi_3];

% create a vector of y
yv = [];
for i = 1:size(y,2)
    yv = [yv; y(:,i)];
end

theta = pinv(phi'*phi)*phi'*yv;

% [Ux,Sx,Vx] = svd(phi);
% S1 = Sx(1:8,1:8);
% V1 = Vx(1:10,1:8);
% U1 = Ux(1:1300,1:8);
% theta = V1/S1*U1'*yv;

xt0 = theta(1:n);
Bt = reshape(theta(n+1:3*n),[n 2]);
Dt = reshape(theta(3*n+1:end),[Ly,Lu]);

toc
%% Create the system
sys = ss(At,Bt,Ct,Dt);
G = tf(sys);
xs = zeros(n,numel(time));

xs(:,1) = xt0;
for i = 1:numel(time)-1
   xs(:,i+1) = At*xs(:,i) + Bt*u(:,i);
end

ys = Ct*xs+Dt*u;

%% Save system matrices in a struct
SID.At = At;
SID.Bt = Bt;
SID.Ct = Ct;
SID.Dt = Dt;
SID.xt0 = xt0;

save('..\..\Data\SID.mat','SID');
%% Plotting

compSubID = figure('Name','SubID Comparison');

subplot(2,1,1)
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
title('Temperature Heater 1')
hold on
plot(time,y1,'b.')
plot(time,ys(1,:),'r','Linewidth',2)
hold off
legend({'y_{1,exp}','y_{1,subID}'},'Location','northeast')

subplot(2,1,2)
xlabel('Time in s')
ylabel('Temperature in /degC',  'Interpreter', 'Latex')
title('Temperature Heater 2')
hold on
plot(time,y2,'b.')
plot(time,ys(2,:),'r','Linewidth',2)
hold off
legend({'y_{2,exp}','y_{2,subID}'},'Location','northeast')

sgtitle('Comperison experiment data and subspace identification output');

% save figure
compSubID = gcf;
compSubID.Renderer = 'painters';
saveas(compSubID, '..\..\Latex\images\SYSID\compSubID', 'svg');

% Singular values
singularVal = figure('Name','Singular Values Y0sN');

xlabel('singular values')
ylabel('Amplitude')
title('Singular Values of Y0sN')
plot(svd(Y0sN),'ro','LineWidth',1);
set(gca,'YScale','log')

% save figure
singularVal = gcf;
singularVal.Renderer = 'painters';
saveas(singularVal, '..\..\Latex\images\SYSID\singularVal', 'svg');
