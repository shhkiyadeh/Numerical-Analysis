% Sanaz Hami Hassan Kiyadeh, Test2, Question 2:

clc
clear
h=0.001;
n=(1/h)-1;
x=h*(0:n+h);
for i=1:n+1
    if i<=n-1
        Q1(i)=(1/h^2)*integral(@(X)(h.*i+h-X).*(X-h.*i).*exp(X),h.*i,h.*i+h);
        Q2(i)=(1/h)^2*integral(@(X)(X-h.*i+h).^2.*exp(X),h.*i-h,h.*i);
        Q3(i)=(-1/h)^2*integral(@(X)(h.*i+h-X).^2.*exp(X),h.*i,h.*i+h);
        Q4(i)=(1/h)^2*integral(@(X)exp(X),h.*i-h,h.*i);
        Q5(i)=(1/h)*integral(@(X)(X-h.*i+h).*(X+(2-X).*exp(X)),h.*i-h,h.*i);
        Q6(i)=(1/h)*integral(@(X)(h.*i+h-X).*(X+(2-X).*exp(X)),h.*i,h.*i+h);
    elseif i==n
        Q2(i)=(1/h)^2*integral(@(X)(X-h.*i+h).^2.*exp(X),h.*i-h,h.*i);
        Q3(i)=(-1/h)^2*integral(@(X)(h.*i+h-X).^2.*exp(X),h.*i,h.*i+h);
        Q4(i)=(1/h)^2*integral(@(X)exp(X),h.*i-h,h.*i);
        Q5(i)=(1/h)*integral(@(X)(X-h.*i+h).*(X+(2-X).*exp(X)),h.*i-h,h.*i);
        Q6(i)=(1/h)*integral(@(X)(h.*i+h-X).*(X+(2-X).*exp(X)),h.*i,h.*i+h);
    elseif i==n+1
        Q4(i)=(1/h)^2*integral(@(X)exp(X),h.*i-h,h.*i);
    end 
end
for i=1:n
    if i<=n-1
        alpha(i)=Q4(i)+Q4(i+1)+Q2(i)+Q3(i);
        beta(i)=Q1(i)-Q4(i+1);
        b(i)=Q5(i)+Q6(i);
    elseif i==n
        alpha(i)=Q4(i)+Q4(i+1)+Q2(i)+Q3(i);
        b(i)=Q5(i)+Q6(i);
    end 
end
a(1)=alpha(1);
s(1)=beta(1)./alpha(1);
z(1)=b(1)./a(1);
for i=2:n-1
    a(i)=alpha(i)-beta(i-1)*s(i-1);
    s(i)=beta(i)/a(i);
    z(i)=(b(i)-beta(i-1)*z(i-1))/a(i);
end
a(n)=alpha(n)-beta(n-1)*s(n-1);
z(n)=(b(n)-beta(n-1)*z(n-1))/a(n);
c(n)=z(n);
for i=n-1:-1:1
    c(i)=z(i)-s(i)*c(i+1);
end

xi=(h:h:1-h);
for i=1:length(xi) 
    approximate(i)= phix(xi(i),n,h,x,c);
    exact(i)=(xi(i)-1)*(exp(-xi(i))-1);
end
error=sum((approximate-exact).^2); 
e=(approximate-exact).^2;

figure(1)
plot(xi,approximate,'ro','LineWidth',0.5)
hold on
plot(xi,exact,'k-','LineWidth',1)
legend('Approximate','Exact')
xlabel('x')
ylabel('y(x)')
title('Piece-wise Linear Rayleigh-Ritz Method')

figure(2)
plot(xi,e)
xlabel('x')
ylabel('error')
title('Piece-wise Linear Rayleigh-Ritz Method')

