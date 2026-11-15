%BPM_focusing_lens.m
%Simulation of Gaussian Beam Focused by a Lens Using BPM
%Project part 1 (Diffractive Propagation Method by a Pseudo-Spectral Approach)
%Paramters used for simulation :
%Ld (light wavelength) =1.03, w0 (waist)=0.05,
%dz(sample distance along z)=0.1, Z(total final distance away from lens) =30,
% f(focal length)= 10
clear
%Gaussian Beam
    N=512; %N : sampling number
    L=2*10^-3; %Display area (window of 2 mm, about 4.5 times the beam width at Z)
    %Ld=input('wavelength of light in [micrometers] = ');
    Ld=1.03;Ld=Ld*10^-6;
    k0=(2*pi)/Ld;
    %w0=input('Waist of Gaussian Beam in [mm] = ');
    w0=0.05;w0=w0*10^-3;
    %dz=input('step size of z (dz) in [mm] = ');
    dz=0.1;dz=dz*10^-3;
    %Z=input('destination of z in [mm] = ');
    Z=30;Z=Z*10^-3;
%Focal length of Lens
    %f=input('Focal length of lens in [mm]= ');
    f=10;f=f*10^-3;
%dx : step size
    dx=L/N;dy=L/N;
    for n= 1:N
      for m=1:N
     %Space axis
        x(m)=(m-1)*dx-L/2;
        y(n)=(n-1)*dy-L/2;
     %Frequency axis
        Kx(m)=(2*pi*(m-1))/(N*dx)-((2*pi*N)/(N*dx))/2;
        Ky(n)=(2*pi*(n-1))/(N*dy)-((2*pi*N)/(N*dy))/2;
      end
    end

[X,Y]=meshgrid(x,y);
[KX,KY]=meshgrid(Kx,Ky);

%Gaussian Beam in space domain (normalisation of the project sheet)
Gau_ini=1/(w0*pi^.5)*exp(-(X.^2+Y.^2)./(w0^2));

%Energy of the initial Gaussian beam
Energy_ini=dx*dy*sum(sum(abs(Gau_ini).^2))

%Lens Equation (named Lens so that the display area L is not overwritten)
Lens=exp(j*k0/(2*f)*(X.^2+Y.^2));

%Gaussian Beam passed through the lens
Gau_L=Gau_ini.*Lens;

%Free space transfer function of propagation
H=exp(j/(2*k0)*dz*(KX.^2+KY.^2));

%Iterative Loop
Gau=Gau_L;
n1=0;
for z=0:dz:Z-dz  %stops at Z-dz so that the last plane is exactly Z
    n1=n1+1 ;
    Zp(n1)=z+dz;
    %Gaussian Beam in Frequency domain
      FGau=fft2(Gau);
    %Propagated Gaussian beam in Frequency domain
      FGau=FGau.*fftshift(H);
    %Propagated Gaussian beam in space domain
      Gau=ifft2(FGau);
    %Step propagation through medium
      Gau_pro(:,n1)=Gau(:,N/2+1);  %central column (x=0)
end

%Energy of the final propagated Gaussian beam,
%to check conservation of energy
Energy_pro=dx*dy*sum(sum(abs(Gau).^2))
Energy_error=abs(Energy_pro-Energy_ini)/Energy_ini

%Beam width along z (1/e of the maximum amplitude, linear interpolation)
Gau_pro=[Gau_L(:,N/2+1) Gau_pro];  %add the lens plane z=0
Zp=[0 Zp];
for l=1:length(Zp)
    P=abs(Gau_pro(:,l))';
    [Pmax,im]=max(P);
    iL=find(P(1:im)<Pmax/exp(1),1,'last');
    iR=im-1+find(P(im:end)<Pmax/exp(1),1,'first');
    yL=interp1(P(iL:iL+1),y(iL:iL+1),Pmax/exp(1));
    yR=interp1(P(iR-1:iR),y(iR-1:iR),Pmax/exp(1));
    W(l)=(yR-yL)/2;
end

%Gaussian beam theory (waist at the lens)
zR=pi*w0^2/Ld;              %Rayleigh length before the lens
d_f=f/(1+(f/zR)^2);         %position of the focus
w_f=w0/sqrt(1+(zR/f)^2);    %waist at the focus
zR_f=pi*w_f^2/Ld;           %Rayleigh length after the lens
W_th=w_f*sqrt(1+((Zp-d_f)/zR_f).^2);

%Comparison BPM / theory
[W_min,i_min]=min(W);
Focus_BPM_mm=Zp(i_min)*10^3
Focus_theory_mm=d_f*10^3
Waist_BPM_um=W_min*10^6
Waist_theory_um=w_f*10^6
Width_Z_BPM_um=W(end)*10^6
Width_Z_theory_um=W_th(end)*10^6
Max_relative_deviation=max(abs(W-W_th)./W_th)

%axis in mm scale
x=x*10^3;
y=y*10^3;
Zp=Zp*10^3;
MAX1=max(max(abs(Gau_ini)));
MAX2=max(max(abs(Gau)));
MAX=max([MAX1 MAX2]);
Zoom=0.5; %half width of the displayed area in mm

figure(1);
mesh(x,y,abs(Gau_ini))
title('Initial Gaussian Beam')
xlabel('x [mm]')
ylabel('y [mm]')
axis([-Zoom Zoom -Zoom Zoom 0 MAX])
axis square

figure(2);
mesh(x,y,abs(Gau))
title('Propagated Gaussian Beam at Z')
xlabel('x [mm]')
ylabel('y [mm]')
axis([-Zoom Zoom -Zoom Zoom 0 MAX])
axis square

figure(3)
for l=1:10:length(Zp)  %one profile every 1 mm
  plot3(y',Zp(l)*ones(size(y')),abs(Gau_pro(:,l)))
  hold on
end
axis([-Zoom Zoom min(Zp) max(Zp)])
grid on
title('Beam profile along z')
xlabel('x [mm]')
ylabel('z [mm]')
hold off

figure(4)
A=max(abs(Gau_pro));
B=diag(1./A);
N_Gau_pro=abs(Gau_pro)*B;
contour(Zp,y,N_Gau_pro, [ exp(-1) exp(-1)], 'k')
axis([min(Zp) max(Zp) -Zoom Zoom])
grid on
title('Beam waist along z')
xlabel('z [mm]')
ylabel('x [mm]')

figure(5)
subplot(1,2,1)
imagesc(x,y,abs(Gau_L).^2)
axis([-Zoom Zoom -Zoom Zoom])
axis square
colorbar
title('Intensity at the lens plane (z=0)')
xlabel('x [mm]')
ylabel('y [mm]')
subplot(1,2,2)
imagesc(x,y,abs(Gau).^2)
axis([-Zoom Zoom -Zoom Zoom])
axis square
colorbar
title('Intensity at Z')
xlabel('x [mm]')
ylabel('y [mm]')

figure(6)
plot(Zp(1:5:end),W(1:5:end)*10^6,'bo',Zp,W_th*10^6,'k-',Zp,w0*abs(1-Zp/(f*10^3))*10^6,'r--')
grid on
legend('BPM (1/e of max)','Gaussian beam theory','Geometrical optics','Location','northwest')
title('Beam width along z')
xlabel('z [mm]')
ylabel('w(z) [\mum]')
%Beam profile before propagation (z=0) : before and just after the lens
Initial_width_um=W(1)*10^6  %1/e half-width at z=0 (should be w0 = 50 um)
P_ini=abs(Gau_ini(:,N/2+1));          %amplitude before the lens (central cut)
P_L=abs(Gau_L(:,N/2+1));              %amplitude just after the lens
Phi_L=unwrap(angle(Gau_L(:,N/2+1)));  %phase just after the lens
Phi_L=Phi_L-Phi_L(N/2+1);             %phase referenced to the centre

figure(7)
subplot(1,3,1)  %intensity map at z=0
imagesc(x,y,abs(Gau_ini).^2)
axis([-Zoom/2 Zoom/2 -Zoom/2 Zoom/2])
axis square
colorbar
title('|E|^2 before propagation (z=0)')
xlabel('x [mm]')
ylabel('y [mm]')
subplot(1,3,2)  %amplitude profile : unchanged by the lens
plot(y,P_ini,'b',y,P_L,'r--')
hold on
plot([-w0 w0]*10^3,max(P_ini)/exp(1)*[1 1],'k-o')  %1/e width = w0
hold off
axis([-Zoom/2 Zoom/2 0 1.1*max(P_ini)])
grid on
legend('before the lens','after the lens','1/e width = w_0')
title('Amplitude profile at z=0 (x=0)')
xlabel('y [mm]')
ylabel('|E|')
subplot(1,3,3)  %phase profile : flat before, parabolic after the lens
plot(y,zeros(size(y)),'b',y,Phi_L,'r',y,k0*(y*10^-3).^2/(2*f),'k:')
axis([-Zoom/2 Zoom/2 -1 1.1*k0*(Zoom/2*10^-3)^2/(2*f)])
grid on
legend('before the lens (flat)','after the lens','k_0 y^2 / 2f')
title('Phase profile at z=0 (x=0)')
xlabel('y [mm]')
ylabel('phase [rad]')
