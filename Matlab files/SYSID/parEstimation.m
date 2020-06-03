%% StateSpace system identification, using parameter fitting
% Add some paths
addpath('..\')
addpath('..\')

% Include tclab.m for initialization
tclab;

% Load date
load('Data\blockResponse.mat');
%time.block(end) = [];
%Input corresponding to the data
u = zeros(2,length(time.block));
u(1,1:50) = 50;
u(2,1:50) = 0;
u(1,51:120) = 25;
u(2,51:120) = 50;
u(1,121:200) = 10;
u(2,121:200) = 40;
u(:,201:260) = 40;
u(:,261:350) = 0;
u(1,351:400) = 25;
u(2,351:400) = 50;
u(1,401:end) = 50;
u(2,401:end) = 25;
%%
emissivity = 0.9;
sigma = 5.67e-8;
Tinf = 23+273;
Cp = 500;
mass = 4e-3;
As = 2e-4;
A = 1e-3;



%States and their derivative
syms TH1 TH2 TC1 TC2 'real'
%syms TH1d TH2d TC1d TC2d 'real'

%Inputs
syms Q1 Q2 'real'

%Parameters to estimate
syms U  Us alpha1 alpha2 tau 'real'
theta = [U; Us; alpha1; alpha2; tau];
QC12 = Us*As*(TH2-TH1);
QR12 = emissivity*sigma*(TH2^4 - TH1^4);

%Dynamics
TH1d = (U*A*(Tinf-TH1) + emissivity*sigma*A*(Tinf^4 - TH1^4) + QC12 + QR12 + alpha1*Q1)/(mass*Cp);
TH2d = (U*A*(Tinf-TH2) + emissivity*sigma*A*(Tinf^4 - TH2^4) - QC12 - QR12 + alpha2*Q2)/(mass*Cp);
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
C = [0, 0, 1, 0; 0, 0, 0, 1]; %C matrix is known; state 3 and 4 are measured
D = zeros(2,2);

%Initial parameters geuss, and linearise around 30C:
% p = 20; %number of simulations that will be performed
% U_ = linspace(1,20,p); %1 is the lower bound, 20 is the upperbound
% Us_ = linspace(5,40,p);


theta_ = [10, 20, 1, 0.75, 10]';
A_ = double(subs(A,[theta; TH1; TH2],[theta_; 273+30; 273+30]));
B_ = double(subs(B,theta,theta_));

%Create the state space model
linsys = ss(A_,B_,C,D);
%Simulate
time.sim = linspace(0,time.block(end),length(time.block));
t1.inter = interp1(time.block,t1.block,time.sim);
t2.inter = interp1(time.block,t2.block,time.sim);
[yhat,~,xsim] = lsim(linsys,u,time.sim);

% %Cost function
% %Set up phi for very timestep:
% A_power(1).A = eye(size(A));
% for i = 1 : length(time.block)
%     sum = zeros(2,4);
%     if i > 1
%         for j = 0: i-1
%             sum = sum + kron(u(:,i)',C*A_power(i-1-j).A);
%         end
%     end
%     phi(i).mat = [C*A_power(i).A, sum, kron(u(:,i)',eye(2))]; 
%     A_power(i+1).A = A*A_power(i).A;
% end