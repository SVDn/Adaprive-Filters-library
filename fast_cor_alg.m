function [Zout,Hout]=fast_cor_alg(Sin, Des_sig, mu, N, M)
% Реализация быстрого корреляционного алгоритма
% Версия 0.1
% Реализация для дейсвительных сигналов
% Входные данные: Входной сигнал, сигнал отклика системы
% Параметры: mu,N,M
% Выходные данные: остаточный эхосигнал,импульсная характеристика
% ПРОВЕРИТЬ ВСЕ ИНДЕКСЫ!!! Есть сомнения, что все правильно.
% Для сдвига массивов используется пересылка из массива в массив, что не есть хорошо. Реализовать циклический буфер.  


B=0.1;


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


%Переменные для быстрого вычисления дисперсии и функции ВКФ
Dx=zeros(1,M); % хранение отсчетов дисперсии входноо сигнала 
Dx_sup=zeros(1,M); % вспомогательный буфер (это не есть правильно, потом переделать)
Disp=B;

Pz=zeros(N,M); % хранение отсчетов ВКФ
Pz_sup=zeros(1,M); % вспомогательный буфер (это не есть правильно, потом переделать)
VKF=zeros(1,N);



for k = 1:Len
   
    % быстрое вычисление дисперсии (с самого начала, чтобы к моменту N+N уже была норм. оценка)
   for j = 2:M
         Dx_sup(j)=Dx(j-1);       
   end
   for j = 2:M
          Dx(j)=Dx_sup(j);       
   end
   Dx(1)= x(k)*x(k); 
   Disp = Disp + Dx(1)-Dx(M);   
     
        
      
    if k<=N+N % N+N время накопления начальных данных
        z(k)=x(k);
    else
  % Расчте оклика фильтра
     s(k) = 0;
     for i =1:N
         s(k)= s(k) + x(k-i)*ha(i);
     end
         
     z(k)=d(k)- s(k);
         
     
  %Расчет ВКФ быстрым способом (с N+N, чтобы к моменту N+N уже было достатчное количество отсчетов)
     for i = 1:N   
         for j = 2:M
               Pz_sup(j)=Pz(i,j-1);  
         end
         for j = 2:M              
               Pz(i,j)=Pz_sup(j);
         end
          Pz(i,1)= x(k-i)*z(k); 
          VKF(i) = VKF(i) + Pz(i,1)-Pz(i,M); 
     end 
     
     
     
     
     %Расчет испульсной характеристики фильтра
     for i = 1:N       
         ha (i) = ha(i) + (mu*VKF(i))/Disp; 
     end
    end
        
end

Zout=z;
Hout=ha;


end

