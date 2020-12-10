
% Пример использования метода корреляции виртуальных сигналов

%clear all
% Initialization
mu =0.0000005;            % MCVS step size
N=32;               % filter length
M=32;               % correlation function length
B=0.001;


% input signal
Len=400000;           % input signal length
% x=wgn(1,Len,0);     % input signal
% noise=wgn(1,Len,0);

% here is QAM64 signals
 x=transpose(QAM64_q(1:Len));
 noise=transpose(QAM64_i(1:Len));


% massive initialisation
%ha = fir1(31,0.5)*0; % filter init

% FIR system to be identified
%b  = fir1(31,0.5);     % FIR system to be identified

% FIR system to be identified
b=[0,0,0,0,0,0.1,0,0,0.9,0.1,0.1,0,0.1,0,0,0,0,0,0,0];

% Desired signal or echo signal
d0  = filter(b,1,x);  % Desired signal
d = d0 + noise; % Desired signal plus noise


[z,H]=improved_mcvs_for_CnC(x, d, mu, N, M);


zwn = z - noise(1:Len);



v=1;
for v= 1:10
%calculate of the signal's level
y1=0;
y4=0;
for k = drange(Len-2000*v:Len-2000*(v-1))
y1 = y1+ zwn(k)^2;
y4 = y4+noise(k)^2;
end
y1 = y1/2000;
y4 = y4/2000;
L1 = 10*log10(y4/y1)
end

% Plot the results
subplot(1,1,1);
%plot (1:Len, x(1:Len),'-.', 1:Len, noise(1:Len),':', 1:Len, zwn,'-');
plot (1:Len,zwn,'-');
legend('Input', 'Noise', 'Error'); % plot without noise
title('System Identification of FIR Filter');
xlabel('Time Index'); ylabel('Signal Value');



% %Plot the results
% subplot(2,1,1);
% plot(1:Len,[x; z]); % plot without noise
% legend('Input', 'Error'); % plot without noise
% title('System Identification of FIR Filter');
% xlabel('Time Index'); ylabel('Signal Value');


% % Pilse response
% subplot(2,1,2);
% stem(1:N,ha);
% legend('Actual','Estimated'); grid on;
% xlabel('Coefficient #'); ylabel('Coefficient Value');


