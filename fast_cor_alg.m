function [Zout,Hout]=fast_cor_alg(Sin, Des_sig, mu, N, M)
% ���������� �������� ��������������� ���������
% ������ 0.1
% ���������� ��� ������������� ��������
% ������� ������: ������� ������, ������ ������� �������
% ���������: mu,N,M
% �������� ������: ���������� ���������,���������� ��������������
% ��������� ��� �������!!! ���� ��������, ��� ��� ���������.
% ��� ������ �������� ������������ ��������� �� ������� � ������, ��� �� ���� ������. ����������� ����������� �����.  


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


%���������� ��� �������� ���������� ��������� � ������� ���
Dx=zeros(1,M); % �������� �������� ��������� ������� ������� 
Dx_sup=zeros(1,M); % ��������������� ����� (��� �� ���� ���������, ����� ����������)
Disp=B;

Pz=zeros(N,M); % �������� �������� ���
Pz_sup=zeros(1,M); % ��������������� ����� (��� �� ���� ���������, ����� ����������)
VKF=zeros(1,N);



for k = 1:Len
   
    % ������� ���������� ��������� (� ������ ������, ����� � ������� N+N ��� ���� ����. ������)
   for j = 2:M
         Dx_sup(j)=Dx(j-1);       
   end
   for j = 2:M
          Dx(j)=Dx_sup(j);       
   end
   Dx(1)= x(k)*x(k); 
   Disp = Disp + Dx(1)-Dx(M);   
     
        
      
    if k<=N+N % N+N ����� ���������� ��������� ������
        z(k)=x(k);
    else
  % ������ ������ �������
     s(k) = 0;
     for i =1:N
         s(k)= s(k) + x(k-i)*ha(i);
     end
         
     z(k)=d(k)- s(k);
         
     
  %������ ��� ������� �������� (� N+N, ����� � ������� N+N ��� ���� ���������� ���������� ��������)
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
     
     
     
     
     %������ ���������� �������������� �������
     for i = 1:N       
         ha (i) = ha(i) + (mu*VKF(i))/Disp; 
     end
    end
        
end

Zout=z;
Hout=ha;


end

