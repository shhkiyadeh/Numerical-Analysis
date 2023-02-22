function ph = phix(v,n,h,x,c)
for i=1:n
    if     v>=0     && v<=x(i)  ; f(i)=0;
    elseif v>x(i)   && v<=x(i+1); f(i)=(v-x(i))/h;
    elseif v>x(i+1) && v<=x(i+2); f(i)=(x(i+2)-v)/h;
    elseif v>x(i+2) && v<=1;      f(i)=0;
    end
end
ph=c*f';
end
