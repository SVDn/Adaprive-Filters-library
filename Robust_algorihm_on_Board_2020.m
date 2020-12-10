
function [Zout,Hout]=Robust_algorihm_on_Board_2020(Sin, Des_sig, mu, N, M, lam)
% Робастный алгоритм, усоврешенствованная версия, представленная на
% конференции on Board 2020 (МТУСИ)



% input signal
Len=length(Sin);           % input signal length

% here is QAM64 signals
x=transpose(Sin(1:Len));
d=transpose(Des_sig(1:Len));

%initialisation of massives
s=zeros(1,Len);
z=zeros(1,Len);
s_1=zeros(1,Len);
Pz=zeros(1,Len);
z_1=zeros(1,Len);
zwn=zeros(1,Len);
ha = zeros(1,N); % filter init
B=0.1;

% Pulse response


% [X_mul,R] = corrmtx(x(end-1000:end),N-1,'covariance');
% inR=inv(R);

[X_mul,R] = corrmtx(x(1:10000),N-1,'covariance');
inR=inv(R);

for k = 2*N+M+1:Len
    
%     if mod (k,10000)==0
%       [X_mul,R] = corrmtx(x(k-9999:k),N-1,'covariance');
%        inR=inv(R);
%     end
    
    
% Расчет вектора Калмана и обратной АКМ
%             XL = x(k-N+1:k);
%             inRXL = inR*XL;
%             invden = 1/(lam + XL'*inR*XL);
%             g = inRXL*invden;              %  вектор Калмана
%             inR = (inR - g*XL'*inR)/lam;   %  обратная автокорреляционная матрица

     % calculation of signal filter
     s(k) = 0;
     for i =1:N
         s(k)= s(k) + x(k-i+1)*ha(i);
     end
         
     z(k)=d(k)- s(k);
    
       % Использование обычного остаточного эхосигнала 
     % VKF and coefficients   
      for i = 1:N
         VKF=0;
         Dx=0;
         Dz=0;
         for j = 1:M
             VKF=VKF+x(k-j+1-i+1)*z(k-j+1);
             Dx=Dx+x(k-j+1-i+1)^2; % можно убрать, заменив на мю 
             Dz=Dz+z(k-j+1)^2; % можно убрать, заменив на мю
         end
         %Pz (i) = VKF/(sqrt(Dx*Dz)); 
         Pz (i) = VKF/M;
      end

   % Использование виртуального остаточного эхосигнала   
%  % virtual signals for adfptation
%      for j=1:M+N
%          s_1(j)=0;
%          for i =1:N
%          s_1(j)= s_1(j) + x(k-i+1-j+1)*ha(i);
%          end
%          z_1(j)=d(k-j+1)- s_1(j);
%      end
%        
%      % VKF and coefficients   
%       for i = 1:N
%          VKF=0;
%          Dx=0;
%          Dz=0.01;
%          for j = 1:M
%              VKF=VKF+x(k-j+1-i+1)*z_1(j);
%              Dx=Dx+x(k-j+1-i+1)^2; % можно убрать, заменив на мю 
%              Dz=Dz+z_1(k-j+1)^2; % можно убрать, заменив на мю
%          end
%          %Pz (i) = VKF/(sqrt(Dx*Dz));
%          Pz (i) = VKF;
%       end        


      deltaH=inR*transpose(Pz(1:N));
      
       for i = 1:N
         ha (i) = ha(i) + mu*deltaH(i)*(abs(transpose(Pz(i)))^2); 
       end
        
end

Zout=z;
Hout=ha;

end


