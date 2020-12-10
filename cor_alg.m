function [Zout,Hout]=cor_alg(Sin, Des_sig, mu, N, M)
% Реализация корреляционного алгоритма
% Версия 0.1
% Реализация для дейсвительных сигналов
% Входные данные: Входной сигнал, сигнал отклика системы
% Параметры: mu,N,M
% Выходные данные: остаточный эхосигнал,импульсная характеристика
% ПРОВЕРИТЬ ВСЕ ИНДЕКСЫ!!! Есть сомнения, что все правильно.


B=0.00001;


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
    
    if k<=N+M
        z(k)=x(k);
    else
  % calculation of signal from filter
     s(k) = 0;
     for i =1:N
         s(k)= s(k) + x(k-i)*ha(i);
     end
         
     z(k)=d(k)- s(k);
    
    Dx=B;
    for i = 1:M
       Dx = Dx + x(k-i+1)*x(k-i+1);  
    end
     %Dx = Dx/M;
        
     
     for i = 1:N
         VKF=0;
         for j = 0:M
             VKF=VKF+x(k-j-i)*z(k-j);
         end
         ha (i) = ha(i) + (mu*VKF)/Dx; 
     end
    end
        
end

Zout=z;
Hout=ha;


end

