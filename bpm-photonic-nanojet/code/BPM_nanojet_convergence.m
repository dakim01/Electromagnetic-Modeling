
%BPM_nanojet_convergence.m
%Numerical checks of the 2D BPM nanojet (silica cylinder in water)
%Project part 2, task 5 : convergence in dx, dz and window size, energy, absorbing boundary
%Paramters of the reference simulation :
%Ld =1.03, D =10.3, w0 =39.8, np =1.450, nbg =1.326, N=4096, L=160, dz=0.025 (micrometers)
%Each line of Runs is one simulation : [N  L(um)  dz(um)  absorber]
%absorber : 0 = none, 1 = super-Gaussian mask applied at every step (first version),
%           2 = graded absorber defined per unit length (independent of dz)
clear
%Gaussian Beam
    Ld=1.03;Ld=Ld*10^-6;
    k0=(2*pi)/Ld;
    w0=39.8;w0=w0*10^-6;
%Microcylinder
    D=10.3;D=D*10^-6;R=D/2;
    zc=5*10^-6+R;
    Zw=zc;
    Z=zc+R+40*10^-6;
    np=1.450;nbg=1.326;
    kbg=k0*nbg;
%Absorbing boundary (graded) : starts at 80% of the half window
    alpha=5*10^5;  %absorption coefficient at the window edge [1/m]
%List of simulations
    Runs=[1024 160 0.025  2;   %--- transverse step dx
          2048 160 0.025  2;
          4096 160 0.025  2;
          8192 160 0.025  2;
          4096 160 0.2    2;   %--- longitudinal step dz
          4096 160 0.1    2;
          4096 160 0.05   2;
          4096 160 0.0125 2;
          3072 120 0.025  2;   %--- window size (dx fixed = 39 nm)
          6144 240 0.025  2;
          4096 160 0.025  0;   %--- absorbing boundary
          4096 160 0.025  1];
    names={'dx','dx','dx','dx','dz','dz','dz','dz','L','L','no absorber','mask every step'};

for r=1:size(Runs,1)
    N=Runs(r,1);L=Runs(r,2)*10^-6;dz=Runs(r,3)*10^-6;abs_type=Runs(r,4);
    %dx : step size
    dx=L/N;
    clear x Kx Zp I_axis P_pro
    for m=1:N
     %Space axis
        x(m)=(m-1)*dx-L/2;
     %Frequency axis
        Kx(m)=(2*pi*(m-1))/(N*dx)-((2*pi*N)/(N*dx))/2;
    end
    ic=N/2+1;
    ix=find(abs(x)<=4*10^-6);
    %Absorbing boundary
    if abs_type==0
        Abs=ones(1,N);
    elseif abs_type==1
        Abs=exp(-(abs(x)/(0.42*L)).^30);
    else
        xa=0.8*L/2;
        Abs=exp(-alpha*dz*(max(abs(x)-xa,0)/(L/2-xa)).^2);
    end
    %Gaussian Beam in space domain at z=0
    Gau_w=exp(-x.^2./(w0^2));
    Gau_ini=ifft(fft(Gau_w).*fftshift(exp(j/(2*kbg)*(-Zw)*Kx.^2)));
    Energy_ini=dx*sum(abs(Gau_ini).^2);
    %Free space transfer function of propagation
    H=exp(j/(2*kbg)*dz*Kx.^2);
    %Iterative Loop
    Gau=Gau_ini;
    n1=0;
    Energy_abs=0;
    Nz=round(Z/dz);
    for n=1:Nz
        z=(n-1)*dz;
        n1=n1+1 ;
        Zp(n1)=z+dz;
        %Gaussian Beam in Frequency domain
          FGau=fft(Gau);
        %Propagated Gaussian beam in Frequency domain
          FGau=FGau.*fftshift(H);
        %Propagated Gaussian beam in space domain
          Gau=ifft(FGau);
        %Refractive index map of the slice z+dz/2
          dn=(np-nbg)*((x.^2+(z+dz/2-zc)^2)<=R^2);
        %Step propagation through medium
          Gau=Gau.*exp(-j*k0*dn*dz);
        %Absorbing boundary
          E_before=dx*sum(abs(Gau).^2);
          Gau=Gau.*Abs;
          Energy_abs=Energy_abs+E_before-dx*sum(abs(Gau).^2);
        %Stored intensity
          I_axis(n1)=abs(Gau(ic))^2;
          P_pro(:,n1)=abs(Gau(ix)').^2;
    end
    Energy_pro=dx*sum(abs(Gau).^2);

    %Figures of merit
    after=find(Zp>zc+R);
    [Imax,k]=max(I_axis(after));kmax=after(k);
    %sub-step position of the maximum (parabola through 3 points)
    y0=I_axis(kmax-1);y1=I_axis(kmax);y2=I_axis(kmax+1);
    off=0.5*(y0-y2)/(y0-2*y1+y2);
    d_jet=(Zp(kmax)+off*dz-(zc+R))*10^6;
    Imax=y1-0.25*(y0-y2)*off;
    %FWHM at the maximum
    xs=x(ix);P=P_pro(:,kmax)';i0=find(xs==0);
    iL=find(P(1:i0)<P(i0)/2,1,'last');
    iR=i0-1+find(P(i0:end)<P(i0)/2,1,'first');
    xL=interp1(P(iL:iL+1),xs(iL:iL+1),P(i0)/2);
    xR=interp1(P(iR-1:iR),xs(iR-1:iR),P(i0)/2);
    %length above Imax/2
    k1=kmax;k2=kmax;
    while k1>after(1) && I_axis(k1-1)>=Imax/2, k1=k1-1; end
    while k2<length(Zp) && I_axis(k2+1)>=Imax/2, k2=k2+1; end
    %field left at the window edge (wrap-around indicator)
    Edge=max(abs(Gau([1:round(N/40) N-round(N/40)+1:N])))/max(abs(Gau));

    Res(r,:)=[dx*10^9 dz*10^9 L*10^6 d_jet Imax (xR-xL)*10^6 (Zp(k2)-Zp(k1))*10^6 ...
              100*Energy_abs/Energy_ini abs(Energy_pro+Energy_abs-Energy_ini)/Energy_ini Edge];
    Axis_r{r}=[(Zp-(zc+R))*10^6;I_axis];
    fprintf('run %2d (%-15s) dx=%5.1f nm dz=%5.1f nm L=%3.0f um | d=%6.3f um eta=%6.3f FWHM=%5.3f um L1/2=%6.2f um | absorbed=%7.4f %% balance=%8.1e edge=%8.1e\n', ...
            r,names{r},Res(r,1),Res(r,2),Res(r,3),Res(r,4),Res(r,5),Res(r,6),Res(r,7),Res(r,8),Res(r,9),Res(r,10));
end

%Relative change with respect to the finest run of each family [%]
Ref_dx=Res(4,:);Ref_dz=Res(8,:);Ref_L=Res(10,:);
Change_dx=100*(Res(1:4,4:7)-Ref_dx(4:7))./Ref_dx(4:7)  %columns: d, eta, FWHM, L1/2
Change_dz=100*(Res([5 6 7 3 8],4:7)-Ref_dz(4:7))./Ref_dz(4:7)
Change_L=100*(Res([9 3 10],4:7)-Ref_L(4:7))./Ref_L(4:7)

figure(1)
semilogx(Res(1:4,1),Change_dx,'o-')
grid on
legend('d','\eta','FWHM','L_{1/2}')
title('Convergence with the transverse step dx (L=160 \mum, dz=25 nm)')
xlabel('dx [nm]')
ylabel('change w.r.t. dx=19.5 nm [%]')

figure(2)
semilogx(Res([5 6 7 3 8],2),Change_dz,'o-')
grid on
legend('d','\eta','FWHM','L_{1/2}')
title('Convergence with the propagation step dz (dx=39 nm)')
xlabel('dz [nm]')
ylabel('change w.r.t. dz=12.5 nm [%]')

figure(3)
subplot(1,2,1)
plot(Res([9 3 10],3),Change_L,'o-')
grid on
legend('d','\eta','FWHM','L_{1/2}')
title('Window size (dx=39 nm)')
xlabel('L [\mum]')
ylabel('change w.r.t. L=240 \mum [%]')
subplot(1,2,2)
bar(Res([9 3 10],8))
set(gca,'XTickLabel',{'120','160','240'})
grid on
title('Energy absorbed at the edges')
xlabel('L [\mum]')
ylabel('absorbed energy [%]')

figure(4)
plot(Axis_r{11}(1,:),Axis_r{11}(2,:),'r',Axis_r{12}(1,:),Axis_r{12}(2,:),'g',Axis_r{3}(1,:),Axis_r{3}(2,:),'b--')
axis([-15 40 0 9])
grid on
legend(sprintf('no absorber (edge field %.1e)',Res(11,10)), ...
       sprintf('mask at every step (%.2f %% absorbed)',Res(12,8)), ...
       sprintf('graded absorber (%.2f %% absorbed)',Res(3,8)))
title('Effect of the absorbing boundary on the on-axis intensity')
xlabel('distance from the shadow surface [\mum]')
ylabel('I(0,z) / I_0')
