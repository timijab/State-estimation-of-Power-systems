function [PS] = PS_sparse(H)
[m,n]=size(H);
Ohmega=H*H';
for k=1:m
    for i=1:m
        for j=1:m
            if j~=i
                x(j)=abs(Ohmega(i,k)+Ohmega(j,k));
            end
        end
       mask=x~=0; % mark the non-zero elements
       x=x(mask); % keep the non-zero elements
       x0=sort(x,'ascend');  % order the x so that lomed can be calculated
       y(i)=x0(floor((length(x0)+1)/2)); % calculate the low median 
%         x(x==0)=NaN;
%         y(i)=nanmedian(x);
    end
    y0=sort(y,'ascend'); % order the x so that lomed can be calculated
    sm(k)=1.1926*y0(floor((length(y0)+1)/2));% calculate the low median
%    sm(k)=1.1926*correction_factor(m)*y0(floor((length(y0)+1)/2));%
%    calculate the low median including the correction factor cm
%     y(y==0)=NaN;
%     sm(k)=1.1926*nanmedian(y);
end
for k=1:m
    for i=1:m
        aux(i)=abs(Ohmega(k,i))/sm(i);
    end
    PS(k,1)=max(aux);
end