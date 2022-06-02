% Solve the problem:
n = 120; r = 0.01; % number and radius of disks
c = exp(2i*pi*(1:n)/n); % centers of the disks
rr = r*ones(size(c)); % vector of radii
N = max(0,round(4+.5*log10(r))); % number of terms in expansions
npts = 3*N+2; % number of sample points on circles
circ = exp((1:npts)'*2i*pi/npts); % roots of unity for collocation
z = []; for j = 1:n
z = [z; c(j)+rr(j)*circ]; end % collocation points
A = [0; -ones(size(z))]; % the constant term
zs = 2; % location of the singularity
rhs = [0; -log(abs(z-zs))]; % right-hand side
for j = 1:n
A = [A [1; log(abs(z-c(j)))]]; % the logarithmic terms
for k = 1:N
zck = (z-c(j)).^(-k);
A = [A [0;real(zck)] [0;imag(zck)]]; % the algebraic terms
end
end
X = A\rhs; % solve least-squares problem
e = X(1); X(1) =[]; % constant potential on wires
d = X(1:2*N+1:end); X(1:2*N+1:end) = []; % coeffs of log terms
a = X(1:2:end); b = X(2:2:end); % coeffs of algebraic terms
% Plot the solution:
x = linspace(-1.4,2.2,120); y = linspace(-1.8,1.8,120);
[xx,yy] = meshgrid(x,y); zz = xx+1i*yy; uu = log(abs(zz-zs));
for j = 1:n
uu = uu+d(j)*log(abs(zz-c(j)));
for k = 1:N, zck = (zz-c(j)).^(-k); kk = k+(j-1)*N;
uu = uu+a(kk)*real(zck)+b(kk)*imag(zck); end
end
for j = 1:n, uu(abs(zz-c(j))<rr(j)) = NaN; end
z = exp(pi*1i*(-50:50)'/50);
for j = 1:n, disk = c(j)+rr(j)*z; fill(real(disk),imag(disk),[1 .7 .7])
hold on, plot(disk,'-r'), end
contour(xx,yy,uu,-2:.1:1.2), colormap([0 0 0]), axis([-1.4 2.2 -1.8 1.8])
axis square, plot(real(zs),imag(zs),'.r')