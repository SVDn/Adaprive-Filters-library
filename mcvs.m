
function [Zout,Hout]=mcvs(Sin, Des_sig, mu, N, M)
% Реализация метода корреляции виртаульных сигналов
% Версия 0.1
% Реализация для дейсвительных сигналов
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
s_1=zeros(1,Len);
z_1=zeros(1,Len);
zwn=zeros(1,Len);
ha = zeros(1,N); % filter init


% Pulse response

for k = 2*N+M+1:Len
    

    Dx=B;
    for i = 1:M
       Dx = Dx + x(k-i+1)*x(k-i+1);  
    end
 
    
    
     % calculation of signal filter
     s(k) = 0;
     for i =1:N
         s(k)= s(k) + x(k-i+1)*ha(i);
     end
         
     z(k)=d(k)- s(k);
     
     
     % virtual signals for adfptation
     for j=1:M+N
         s_1(j)=0;
         for i =1:N
         s_1(j)= s_1(j) + x(k-i+1-j+1)*ha(i);
         end
         z_1(j)=d(k-j+1)- s_1(j);
     end
     
     % VKF and coefficients
     for i = 1:N
         VKF=0;
         for j = 1:M
             VKF=VKF+x(k-j+1-i+1)*z_1(j);
         end
         ha (i) = ha(i) + (mu*VKF)/Dx; 
     end
        
end

Zout=z;
Hout=ha;

end


