%% StateSpace system identification, using parameter fitting
% Add some paths
addpath('..\')
addpath('..\')

% Include tclab.m for initialization
%tclab;

% Load date
cd '..\..'
load('Data\blockResponse.mat');

%Input/output corresponding to the data
u = [u1;u2];
y = [t1.block',t2.block'];
%%
emissivity = 0.9;
sigma = 5.67e-8;
Tinf = 23+273.15;
Cp = 500;
mass = 4e-3;
As = 2e-4;
A_par = 1e-3;

T0 = Tinf; %Linearisation Temperature

%States and their derivative
syms TH1 TH2 TC1 TC2 'real'
%syms TH1d TH2d TC1d TC2d 'real'

%Inputs
syms Q1 Q2 'real'

%Parameters to estimate
syms U  Us alpha1 alpha2 tau 'real'
theta = [U; Us; alpha1; alpha2; tau];
QC12 = Us*As*(TH2-TH1);
QR12 = emissivity*sigma*As*(TH2^4 - TH1^4);

%Dynamics
TH1d = (U*A_par*(Tinf-TH1) + emissivity*sigma*A_par*(Tinf^4 - TH1^4) + QC12 + QR12 + alpha1*Q1)/(mass*Cp);
TH2d = (U*A_par*(Tinf-TH2) + emissivity*sigma*A_par*(Tinf^4 - TH2^4) - QC12 - QR12 + alpha2*Q2)/(mass*Cp);
TC1d = (TH1 - TC1)/tau;
TC2d = (TH2 - TC2)/tau;

%State vector and its derivative
x = [TH1; TH2; TC1;TC2];
xd = [TH1d; TH2d; TC1d; TC2d];
%Input vector
input = [Q1; Q2];
%Linearise to get the state space model
A = jacobian(xd, x);
B = jacobian(xd, input);
A = subs(A,[TH1; TH2],[T0; T0]); %Linearisation around T0
C = [0, 0, 1, 0; 0, 0, 0, 1]; %C matrix is known; state 3 and 4 are measured
D = zeros(2,2);

%Initial parameters geuss, and linearise around 30C:
% p = 20; %number of simulations that will be performed
% U_ = linspace(1,20,p); %1 is the lower bound, 20 is the upperbound
% Us_ = linspace(5,40,p);


theta_ = [10, 20, 0.01, 0.005, 10]';
A_ = double(subs(A,theta,theta_));
B_ = double(subs(B,theta,theta_));

%Create the state space model
linsys = ss(A_,B_,C,D);

%Simulate
time.sim = linspace(0,time.block(end),length(time.block));
t1.inter = interp1(time.block,t1.block,time.sim);
t2.inter = interp1(time.block,t2.block,time.sim);
x0 = [0;0;0;0]; %initial state
[yhat,~,xhat] = lsim(linsys,u,time.sim);
xhat = xhat';
%Compute the error vector
EN = zeros(2*length(y),1);
for i = 1:length(y)
    EN(2*i-1:2*i,1) = (y(i,:)-23.15)' - yhat(i,:)';
end
%% Compute Psi(theta)
%To compute Psi(theta), we need to a simulation for every entry of theta
%Compute the derivatives of the state matrices
Psi = zeros(2*length(time.sim),length(theta));
for p = 1:length(theta)
   %Compute the matrix derivatives
   A__ = subs(A,[theta(1:p-1); theta(p+1:end)], [theta_(1:p-1); theta_(p+1:end)]);
   B__ = subs(B,[theta(1:p-1); theta(p+1:end)], [theta_(1:p-1); theta_(p+1:end)]);
   derivative(p).A = double(subs(diff(A__), theta(p), theta_(p)));
   derivative(p).B = double(subs(diff(B__), theta(p), theta_(p)));
   
   %Compute the simulations
   X(p).state(:,1) = x0;
   for i = 1:length(time.sim)
      X(p).state(:,i+1) = A_*X(p).state(:,i) + derivative(p).A*xhat(:,i) + derivative(p).B*u(:,i) ;
      dy(p).state(:,i) = C*X(p).state(:,i);
      
      %Create Psi(theta)
      Psi(2*i-1:2*i,p) = dy(p).state(:,i);
   end
end
Jd = 2/length(time.sim)*Psi'*EN;
H = 2/length(time.sim) * (Psi'*Psi);
theta_est(:,1) = theta_;
for i = 1:10
    theta_est(:,i+1) = theta_est(:,i) - inv(H)*Jd;
end




