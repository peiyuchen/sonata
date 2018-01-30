function [w,accumulate]=dtw(D)
%%
% -----輸入-----
% D：累積成本路徑表
% -----輸出-----
% w：最佳路徑
% accumulate：最佳路徑的累積成本值
%%

n=size(D,1);
m=size(D,2);
N = n;
M = m;
w=[];
w(1,:)=[N,M];
accumulate = D(n,m);
while ((n+m)~=2)
    if (n-1)==0
        m=m-1;
    elseif (m-1)==0
        n=n-1;
    else 
      %[values,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1)]);
      [values,number]=min([D(n-1,m-1),D(n,m-1),D(n-1,m)]);
      switch number
      case 1
        n=n-1;
        m=m-1;
        %accumulate = accumulate + values;
      case 2
        m=m-1;
        %accumulate = accumulate + values;
      case 3
        n=n-1;
        %accumulate = accumulate + values;
      end
    end
    w=cat(1,w,[n,m]);
end
w = flipud(w);