function fn=circonv(x1,x2,N)
% x1=[1 ;2 ;3; 4 ;5];
% x2=[2; 3 ;4];
% N=5;

% x1,x2为输入循环卷积序列，N为循环卷积点数
x2=[x2 ;zeros(N-length(x2),1)];  %补到N点
m=0:N-1;
x=zeros(N,N);
for n=0:N-1
   x(:,n+1)=x2(mod((m-n),N)+1); 
end
Num = length(x1)/N;
for sym =1:Num
fn((sym-1)*N+(1:N))=x*x1((sym-1)*N+(1:N));
end

