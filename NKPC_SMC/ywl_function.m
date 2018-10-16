function[theta]=ywl_learning(y)

N=length(y);
gain=0.05;
% gain=1/N;
num=0;denum=0;
for i=1:N-1
    num=num+gain*(1-gain)^(N-i)*y(i)*y(i+1);
end

for i=1:N
    denum=denum+gain*(1-gain)^(N-i+1)*y(i)*y(i);
end
% 
% for i=1:N-1
%     num=num+gain*(1-gain)^(N-i)*y(i)*y(i+1);
% end


theta=num/denum;
%     theta=sqrt(1-gain)*num/denum;


end