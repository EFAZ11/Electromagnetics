num_point_charges = 200;
x_min = -1;
x_max = 3;
y_min = -2;
y_max = 2;
% distance between equally spaced point charges on unit circle
r = 1 / num_point_charges;
% starting location of exterior point charge
zs_cart = [0; 0];
zs = zs_cart(1) + 1i*zs_cart(2);
n = 1:num_point_charges;
pos = exp(1i*2*pi / num_point_charges * n);
hold on
plot(real(pos), imag(pos), '.r');
plot(real(zs),imag(zs),'.r');
hold off
axis([x_min, x_max, y_min, y_max])
axis square
title('Charged Wire Mesh and Point Charge');
% charges on wires must sum to zero
Aeq = ones(1, num_point_charges);
beq = 0;
H = zeros(num_point_charges, num_point_charges);
f = zeros(num_point_charges, 1);
% arrange point charges in a loop
n = 1:num_point_charges;
pos = exp(1i*2*pi / num_point_charges * n);
% coordinates of external charge
zs_i = zs_cart(1) + 1i*zs_cart(2);
% fill Hessian and gradient with 2-D potentials
% between all charged objects
for row = 1:num_point_charges
    
   H(row,row) = -log(r);    
   f(row) = -log(abs(pos(row) - zs_i));
    
   for col = (row+1):num_point_charges        
       H(row, col) = -log(abs(pos(row) - pos(col)));
       H(col, row) = H(row, col);
   end
end
q_default = quadprog(H, f, [], [], Aeq, beq, [], [], [], ...
            optimoptions(@quadprog, 'Display', 'iter', 'Algorithm', 'interior-point-convex'));
npts = 500;
x = linspace(x_min,x_max,npts);
y = linspace(y_min, y_max,npts);
[X,Y] = meshgrid(x,y); 
Z = X+1i*Y;
zs = zs_cart(1) + 1i*zs_cart(2);
U = log(abs(Z-zs));
n = 1:num_point_charges;
pos = exp(1i*2*pi / num_point_charges * n);
for j = 1:num_point_charges
    U = U + q_default(j)*log(abs(Z-pos(j)));
end
plot(real(pos), imag(pos), '.r');
hold on
contour(X,Y,U,-2:.1:1.2);
axis([x_min, x_max, y_min, y_max]);
axis square
plot(real(zs),imag(zs),'.r');