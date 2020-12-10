
% Пример использования быстрого корреляционного алгоритма

clear all
% Initialization
mu =0.005;            % step size
N=32;               % filter length
M=32;               % correlation function length


% input signal
Len=10000;           % input signal length
x=wgn(1,Len,0);     % input signal

% massive initialisation
ha = fir1(31,0.5)*0; % filter init

% FIR system to be identified
b  = fir1(31,0.5);     % FIR system to be identified

% Desired signal or echo signal
d  = filter(b,1,x);  % Desired signal
%d = d0 + noise; % Desired signal plus noise


[z,H]=fast_cor_alg(x, d, mu, N, M);



% Plot the results
subplot(2,1,1);
plot(1:Len,[x; z]); % plot without noise
legend('Input', 'Error'); % plot without noise
title('System Identification of FIR Filter');
xlabel('Time Index'); ylabel('Signal Value');


% % Pilse response
% subplot(2,1,2);
% stem(1:N,ha);
% legend('Actual','Estimated'); grid on;
% xlabel('Coefficient #'); ylabel('Coefficient Value');


