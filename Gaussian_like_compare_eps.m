% NSDM method applied to a general integral with complex quadratic 
% phase, using both generalized Gauss-Laguerre and Hermite quadrature

close all
omegadata = 2.^(0:7);     % wavelength span
nn = length(omegadata);
d = 10;
f = @(x) exp(-d*x.^2);  % test function
delta = 1;              % distance from the real axis
A = 1+1i; B = -delta; C = 1; % coefficients of the quadratic phase
% A = 1+1i; B = 0; C = 100i;
a = -10; b = 10;        % interval endpoints

fel = zeros(nn,1); I = fel; Iref = fel; I2 = fel; fel2 = fel;   % error initialization
for kk = 1:nn;
    omega = omegadata(kk);
    I(kk) = steepest_descent(omega,A,B,C,a,b,f);    % steepest descent GL
    I2(kk) = steepest_descent_H(omega,A,B,C,a,b,f);    % steepest descent H
    
% compare with matlab's integral: OBSOLETE, we have explicit solution
%     g = @(x) A*x.^2 + B*x + C;
% Iref(kk) = integral(@(x) exp(1i*g(x)*omega).*f(x),...
%     a,b,'RelTol',1e-15,'AbsTol',1e-15);
% Iref(kk) = sqrt(pi)*sqrt(1/omega)*exp(-delta^2*omega/4);
Iref(kk) = sqrt(pi/(d-1i*A*omega))*exp(-B^2/4*omega^2/(d-1i*A*omega) +...
    1i*C*omega);            % exact solution

% x = a:0.001:b;
% plot(x,real(exp(1i*g(x)*omega)),'o-',x,imag(exp(1i*g(x)*omega)),'*-')
% compare to the simpson's rule
% Iref2 = trap(@(x) exp(1i*g(x)*omega).*f(x),a,b,tol);
fel(kk) = abs((Iref(kk)-I(kk))/Iref(kk));
fel2(kk) = abs((Iref(kk)-I2(kk))/Iref(kk));
% disp(fel)
end
epsd = 1./omegadata;
figure
loglog(epsd,abs(I),'o-',epsd,abs(Iref),'*-',...
    epsd,abs(I2),'s-')
legend('|I(GL)|','|I_{ref}|','|I(GH)|','Location','NorthWest')
title(['Integral value as a function of \epsilon, \delta = ',num2str(delta)])
xlabel('\epsilon'); ylabel('I')
% print('-dpdf','Ival_Gausslike_eps.pdf')

figure
loglog(epsd,fel,'o-',epsd, fel2,'*-');%,epsd,epsd.^6)
legend('rel error GL','GH','Location','NorthWest')
title('Relative error between I and I_{ref}')
xlabel('\epsilon'); ylabel('error')
% print('-dpdf','error_Gausslike_eps.pdf')