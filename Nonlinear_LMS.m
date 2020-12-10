function [Zout,Hout]=Nonlinear_LMS(Sin,Sin3,Sin5, Des_sig, mu1, mu3, mu5, N)
% Nonlinear LMS Huawei task



B=0.0001;


% input signal
Len=length(Sin);           % input signal length

% here is QAM64 signals
x=transpose(Sin(1:Len));
x3=transpose(Sin3(1:Len));
x5=transpose(Sin5(1:Len));

d=transpose(Des_sig(1:Len));

%initialisation of massives
s=zeros(1,Len);
s1=zeros(1,Len);
s3=zeros(1,Len);
s5=zeros(1,Len);
z=zeros(1,Len);
zwn=zeros(1,Len);
h1 = zeros(1,N); % linear kernel init
h3 = zeros(1,N); % third kernel init
h5 = zeros(1,N); % 5th kernel init

% Pulse response

for k = 1:Len
    
    if k<=N+N
        z(k)=x(k);
    else
  % calculation of signal from filter
     s(k) = 0;
     s1(k)=0;
     s3(k)=0;
     
     for i =1:N
         s1(k)= s1(k) + x(k-i+1)*h1(i);
     end
     
     for i =1:N
         s3(k)= s3(k) + x3(k-i+1)*h3(i);
     end
     
     for i =1:N
         s5(k)= s5(k) + x5(k-i+1)*h5(i);
     end
     
     
     s(k)=s1(k) + s3(k)+s5(k);
     z(k)=d(k)- s(k);
    
        
    
     for i = 1:N
         mux=x(k-i+1)*z(k);
         h1 (i) = h1(i) + mu1*mux; 
     end
     
     for i = 1:N
         mux=x3(k-i+1)*z(k);
         h3 (i) = h3(i) + mu3*mux; 
     end
     
     for i = 1:N
         mux=x5(k-i+1)*z(k);
         h5 (i) = h5(i) + mu5*mux; 
     end
     
    end
        
end

Zout=z;
Hout=h1; % linear part output 


end

