% File name: Gauss.m
% Propagation of 2D Gaussian pulse by paraxial FT split-step BPM
clear all
close all
N= 10; delta_x = 1/N; minmax=10;
x = -minmax:delta_x:minmax; % creation of space arguments
y = -minmax:delta_x:minmax; % creation of space arguments
delta_z = 5.0; % step size along z-axis
delta_n = 0.4; k_zero = 100;
beta = 20; % propagation constant
%
u_init=(exp(-x.^2))'*exp(-y.^2); % initial Gaussian pulse
mesh(x,y,abs(u_init)); % plots original pulse
 
%
k_x = -minmax:1/N:minmax; % creation of Fourier variables
k_y = -minmax:1/N:minmax; % creation of Fourier variables
%
temp = delta_z/(2*beta);
H_transfer = (exp(1i*temp*k_x.^2))'*(exp(1i*temp*k_y.^2));
H_transfer = fftshift(H_transfer); S_phase = exp(-1i*k_zero*delta_n);
%
for n=1:100
    z = fft2(u_init);
    zz = z.*H_transfer;
    u_prime = ifft2(zz);
    u = S_phase.*u_prime;
    u_init = u;
    mesh(x,y,abs(u_init)) % plots pulse after propagation
    Zlim=1;
    colorbar
    pause(0.1)
end
%
%close all