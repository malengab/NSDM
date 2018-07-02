function I1 = steepest_descent(omega,A,B,C,a,b,f)
g = @(x) A*x.^2 + B*x + C;
% f = @(x) exp(-10*x.^2);%1;
% a = -1; b = 1;
xi = -B/2/A;    % stationary point

% SD paths
uin = @(x,p) sqrt((x-xi).^2 + 1i*p/A);
hx = @(x,p) xi + uin(x,p);
hx2 = @(x,p) xi - uin(x,p);

% derivatives
hprim = @(x,p) 1i/2/A./uin(x,p);
hprim2 = @(x,p) -1i/2/A./uin(x,p);

n = 10;     % nr of points
alpha1 = 0;  % power of generalized G-L
alpha2 = -1/2;
[pj1,w1] = GaussLaguerre(n,alpha1);    % points + weights
pj1 = pj1/omega;
[pj2,w2] = GaussLaguerre(n,alpha2);    % points + weights
pj2 = pj2/omega;

% SD
  F1x = @(x0) -omega^(-1-alpha1).*exp(1i*omega*g(x0))*sum(f(hx(x0,pj1)).*hprim(x0,pj1).*pj1.^(-alpha1).*w1);
  F1xi = @(x0) -omega^(-1-alpha2).*exp(1i*omega*g(x0))*sum(f(hx(x0,pj2)).*hprim(x0,pj2).*pj2.^(-alpha2).*w2);
  F2x = @(x0) -omega^(-1-alpha1).*exp(1i*omega*g(x0))*sum(f(hx2(x0,pj1)).*hprim2(x0,pj1).*pj1.^(-alpha1).*w1);
  F2xi = @(x0) -omega^(-1-alpha2).*exp(1i*omega*g(x0))*sum(f(hx2(x0,pj2)).*hprim2(x0,pj2).*pj2.^(-alpha2).*w2);

I1 = 0;
% pp = 0:0.1:10;

% figure
% select paths
if abs(hx2(a,0)- a)<1e-10;
    I1 = I1 - F2x(a);
%     plot(real(hx2(a,pp)),imag(hx2(a,pp)),'o-')
%     hold all
%         plot(hx(a,pp),'o-')
else
    I1 = I1 - F1x(a);
%     plot(real(hx(a,pp)),imag(hx(a,pp)),'o-')
%     hold all    
%         plot(hx2(a,pp),'o-')
end

if real(hx(xi,0.01)) > real(xi);
    %real(uin(xi,0.001)/delta) > real(xi);
    I1 = I1 + F2xi(xi) - F1xi(xi);
%     plot(real(hx2(xi,pp)),imag(hx2(xi,pp)),'*-')
%     plot(real(hx(xi,pp)),imag(hx(xi,pp)),'s-')
%     
%         plot(hx(xi,pp),'*-')
%     plot(hx2(xi,pp),'s-')
else
    I1 = I1 + F1xi(xi) - F2xi(xi);
%     plot(real(hx(xi,pp)),imag(hx(xi,pp)),'*-')
%     plot(real(hx2(xi,pp)),imag(hx2(xi,pp)),'s-')
%     
%         plot(hx2(xi,pp),'*-')
%     plot(hx(xi,pp),'s-')
end

if abs(hx2(b,0) - b)<1e-10;
    I1 = I1 + F2x(b);
%     plot(real(hx2(b,pp)),imag(hx2(b,pp)),'<-')
%         plot(hx(b,pp),'<-')
else
    I1 = I1 + F1x(b);
%     plot(real(hx(b,pp)),imag(hx(b,pp)),'<-')
%         plot(hx2(b,pp),'<-')
end
% plot(real(xi),imag(xi),'k*',a,0,'k*',b,0,'k*')
% legend('1','2','3','4')
% close
end