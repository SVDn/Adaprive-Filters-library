
% Пример использования быстрого алгоритма корреляции виртуальных сигналов

clear all
% Initialization
mu =0.001;            % MCVS step size
N=32;               % filter length
M=128;               % correlation function length
B=0.001;


% input signal
Len=10000;           % input signal length
x=wgn(1,Len,0);     % input signal


% FIR system to be identified
b  = fir1(31,0.5);     % FIR system to be identified

% Desired signal or echo signal
d  = filter(b,1,x);  % Desired signal
%d = d0 + noise; % Desired signal plus noise


[z,H]=fast_mcvs(x, d, mu, N, M);



% Plot the results
subplot(2,1,1);
plot(1:Len,[x; z]); % plot without noise
legend('Input', 'Error'); % plot without noise
title('System Identification of FIR Filter');
xlabel('Time Index'); ylabel('Signal Value');


% Pilse response
subplot(2,1,2);
plot(1:N,[b; H]);
legend('Actual','Estimated'); grid on;
xlabel('Coefficient #'); ylabel('Coefficient Value');


