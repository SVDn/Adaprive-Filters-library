
% Nonlinear LMS (Huawei task)
% uncomment line 35 to noise

clear all
% Initialization
mu1 =0.003;            % linear kernel step size
mu3= 0.0005;             % third kernel step size
mu5=0.000001;         % 5th kernel step size
N=4;               % filter length
B=0.001;


% input signal
Len=10000;           % input signal length
x=wgn(1,Len,0);     % input signal
noise=wgn(1,Len,-20);     % noise

% system to be identified

h1=[0.1,0.1,0.1];
h3=[0.1,0.1,0.1];
h5=[0.001,0.001,0.001];

x3 = x.*x.*x;
x5 = x.*x.*x.*x.*x;

% Desired signal
d1  = filter(h1,1,x);  % Linear part
d3 = filter(h3,1,x3);  % x^3 part
d5 = filter(h5,1,x5);  % x^5 part

d=d1+d3+d5;

% uncomment line 35 to noise
%d = d + noise; % Desired signal plus noise

[z,H]=Nonlinear_LMS(x,x3,x5, d, mu1, mu3,mu5, N);



% Plot the results
subplot(2,1,1);
plot(1:Len,[x; z]); % plot without noise
legend('Input', 'Error'); % plot without noise
title('Error for nonlinear response');
xlabel('Time Index'); ylabel('Signal Value');


% % Pilse response
% subplot(2,1,2);
% stem(1:N,ha);
% legend('Actual','Estimated'); grid on;
% xlabel('Coefficient #'); ylabel('Coefficient Value');


