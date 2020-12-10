
function [Zout,Hout]=fast_mcvs(Sin, Des_sig, mu, N, M)
% ���������� �������� ��������� ���������� ����������� ��������
% ������ 0.1
% ���������� ��� �������������� ��������
% ������� ������: ������� ������, ������ ������� �������
% ���������: mu,N,M
% �������� ������: ���������� ���������,���������� ��������������
%���� ��� �� �������� - ��������


%�� ������ ������ � �-���� �� 2017 ����

B=0.0001;


% input signal
Len=length(Sin);           % input signal length

% here is QAM64 signals
x=transpose(Sin(1:Len));
d=transpose(Des_sig(1:Len));

%initialization of massives
s=zeros(1,Len);
z=zeros(1,Len);
ha = zeros(1,N); % filter init

%Dx=zeros(1,N)+B;
Pxd=zeros(1,N);
Pxz=zeros(1,N);
R=zeros(N,N);

rxh=zeros(1,N);
delta_ha=zeros(1,N);


% ������������� ������������ ��������
% for k=1:M
%     Dx(j) = Dx(j) + x(i-j)*x(i-j) - x(i-j-M)*x(i-j-M);
%     Pxd(j)=Pxd(j)+x(i-j)*d(i-j)-x(i-j-M)*d(i-j-M);
% end


% ��������� ���������

for i = N+M+2:Len %  ������ k-->i, � ������ ����� ����� ������������ i ��� ����������� ����������� �������
 
    % ������������ ������� �������  
     s(i) = 0;
     for j =1:N
         s(i)= s(i) + x(i-j)*ha(j);
     end         
     z(i)=d(i)- s(i);
    
    % ������ ������� ������� (���� ���� ��� ���� �������)
    Dx=B;
    for j = 1:M
       Dx = Dx + x(i-j)*x(i-j); 
    end
    
    % ������ ���������� �������� ������� � ����������
    for j=1:N
        Pxd(j)=0;
        for k=1:M
            Pxd(j)=Pxd(j)+x(i-j-k-1)*d(i-j); % d(i-j) or d(i-j-k)???
        end
    end
    
    % ������ ������������������ ������� �������� ������� ���� �� ���������
    % M ������ ���� ������ N ��� ���� ������� 
    %[XX,R] = corrmtx(x(i-M:i),N-1
    
    R=zeros(N,N);
    for j=1:M
    XX_add=x(i-N-j+1:i-j)*transpose(x(i-N-j+1:i-j));
    R=R+XX_add;
    end
   
    
    % ��������� ������������������ ������� �� ������� ������ ��
    %R=R;
    rxh=R*transpose(ha);
    
    % ������ ���������� �������� ������� � ����������� ����������
    Pxz=Pxd-mu*transpose(rxh);
    
    % ������ ���������� ������� �������
    
    for j=1:N
        delta_ha(j)=mu*Pxz(j)/Dx;
    end
    % ������������ ������ �������
    ha=ha+delta_ha;  
     
end

Zout=z;
Hout=ha;

end


