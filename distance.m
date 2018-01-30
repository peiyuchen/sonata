% distance - cosine similarity
function di = distance(r1,r2)
%%
% -----��J-----
% r1 & r2�G��V�q
% -----��X-----
% di�G��V�q���Z��(�ۦ���)
%%

r1_zero = 1;
r2_zero = 1;
length = size(r1,1);
for u=1:length
    if r1(u,1) ~= 0
        r1_zero = 0;
    end
    if r2(u,1) ~= 0
        r2_zero = 0;
    end
end
% normalization
%maxr1 = 0;
%maxr2 = 0;
%for j=1:12
    %if r1(j,1) > maxr1
        %maxr1 = r1(j,1);
    %end
    %if r2(j,1) > maxr2
        %maxr2 = r2(j,1);
    %end
%end
%r1 = r1./maxr1;
%r2 = r2./maxr2;
if r1_zero == 1 || r2_zero == 1
    di = 0;
else
    d = dot(r1,r2);
    lengthA = 0;
    lengthB = 0;
    twoA = r1.^2;
    twoB = r2.^2;
    for i=1:length
        lengthA = lengthA + twoA(i,1);
        lengthB = lengthB + twoB(i,1);
    end
    lengthA = sqrt(lengthA);
    lengthB = sqrt(lengthB);
    di=1-(d/(lengthA*lengthB));
end



%dis = (r1-r2).^2;
%ans2 = 0;
%for i=1:12
    %ans2 = ans2 + dis(i,1);
%end
%di = sqrt(ans2);