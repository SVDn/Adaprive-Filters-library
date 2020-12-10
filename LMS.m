function [Zout,Hout]=LMS(Sin, Des_sig, mu, N)
% Входные данные: Входной сигнал, сигнал отклика системы
% Параметры: mu,N,M
% Выходные данные: остаточный эхосигнал,импульсная характеристика
% ПРОВЕРИТЬ ВСЕ ИНДЕКСЫ!!! Есть сомнения, что все правильно.


B=0.0001;


% input signal
Len=length(Sin);           % input signal length

% here is QAM64 signals
x=transpose(Sin(1:Len));
d=transpose(Des_sig(1:Len));

%initialisation of massives
s=zeros(1,Len);
z=zeros(1,Len);
zwn=zeros(1,Len);
ha = zeros(1,N); % filter init

% Pulse response

for k = 1:Len
    
    if k<=N+N
        z(k)=x(k);
    else
  % calculation of signal from filter
     s(k) = 0;
     for i =1:N
         s(k)= s(k) + x(k-i)*ha(i);
     end
         
     z(k)=d(k)- s(k);
    
        
    
     for i = 1:N
         
         mux=x(k-i)*z(k);
         ha (i) = ha(i) + mu*mux; 
     end
    end
        
end

Zout=z;
Hout=ha;


end

