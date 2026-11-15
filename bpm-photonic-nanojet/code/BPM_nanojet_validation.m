
%BPM_nanojet_validation.m
%Validation of the 2D BPM nanojet against the exact Mie solution of a cylinder
%Project part 2, task 4 (validation)
%Exact reference : Mie series for an infinite cylinder under plane-wave illumination,
%E parallel to the cylinder axis (the polarisation closest to the scalar BPM)
%Paramters used for simulation :
%Ld (light wavelength) =1.03, D (cylinder diameter) =10.3, w0 (waist)=39.8,
%dz=0.025, all lengths in micrometers
%4 cases of increasing relative index m=np/nbg :
%silica/water, polystyrene/water, silica/air, n=1.7/air (Chen et al. 2004)
clear
%Gaussian Beam
    N=4096; %N : sampling number
    L=160*10^-6; %Display area
    Ld=1.03;Ld=Ld*10^-6;
    k0=(2*pi)/Ld;  %wavenumber in vacuum, used for the index step
    w0=39.8;w0=w0*10^-6;
    dz=0.025;dz=dz*10^-6;
%Microcylinder
    D=10.3;D=D*10^-6;R=D/2;
    zc=5*10^-6+R;  %centre of the cylinder (front face at z=5 um)
    Zw=zc;         %waist of the beam at the centre of the cylinder
    Z=zc+R+40*10^-6;  %40 um behind the cylinder
%Cases (cylinder index, background index)
    np_list =[1.450 1.572 1.450 1.700];
    nbg_list=[1.326 1.326 1.000 1.000];
    names={'silica / water','polystyrene / water','silica / air','n=1.7 / air (Chen 2004)'};
%dx : step size
    dx=L/N;
    for m=1:N
     %Space axis
        x(m)=(m-1)*dx-L/2;
     %Frequency axis
        Kx(m)=(2*pi*(m-1))/(N*dx)-((2*pi*N)/(N*dx))/2;
    end
ix=find(abs(x)<=4*10^-6);  %stored region |x| < 4 um
xs=x(ix);
ic=find(xs==0);
alpha=5*10^5;xa=0.8*L/2;  %absorbing boundary (graded, per unit length)
Abs=exp(-alpha*dz*(max(abs(x)-xa,0)/(L/2-xa)).^2);

for c=1:length(np_list)
    np=np_list(c);nbg=nbg_list(c);
    kbg=k0*nbg;
    m_rel(c)=np/nbg;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BPM (same algorithm as BPM_nanojet_cylinder.m)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gaussian Beam in space domain at z=0
    Gau_w=exp(-x.^2./(w0^2));
    Gau_ini=ifft(fft(Gau_w).*fftshift(exp(j/(2*kbg)*(-Zw)*Kx.^2)));
    %Free space transfer function of propagation
    H=exp(j/(2*kbg)*dz*Kx.^2);
    %Iterative Loop
    Gau=Gau_ini;
    n1=0;
    clear Zp I_pro
    for z=0:dz:Z-dz
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
          Gau=Gau.*exp(-j*k0*dn*dz).*Abs;
        %Stored intensity
          I_pro(:,n1)=abs(Gau(ix)').^2;
    end
    after=find(Zp>zc+R);
    zB=Zp(after)-(zc+R);       %distance from the shadow surface
    IB=I_pro(ic,after);        %on-axis intensity
    PB_all=I_pro(:,after);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Exact Mie solution (plane wave, E parallel to the axis)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k=kbg;X=k*R;Y=m_rel(c)*X;
    nmax=round(X+4*X^(1/3)+15);
    for n=-nmax:nmax
        Jx=besselj(n,X);Jy=besselj(n,Y);Hx=besselh(n,1,X);
        dJx=(besselj(n-1,X)-besselj(n+1,X))/2;
        dJy=(besselj(n-1,Y)-besselj(n+1,Y))/2;
        dHx=(besselh(n-1,1,X)-besselh(n+1,1,X))/2;
        a(n+nmax+1)=(m_rel(c)*dJy*Jx-Jy*dJx)/(Jy*dHx-m_rel(c)*dJy*Hx);
    end
    %on-axis field behind the cylinder (origin at the centre)
    zM=linspace(R*1.0005,R+40*10^-6,3000);
    EM=exp(j*k*zM);
    for n=-nmax:nmax
        EM=EM+j^n*a(n+nmax+1)*besselh(n,1,k*zM);  %phi=0 on the axis
    end
    IM=abs(EM).^2;
    zM=zM-R;  %distance from the shadow surface

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figures of merit (BPM and Mie)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for s=1:2
        if s==1, zz=zB; II=IB; else, zz=zM; II=IM; end
        [Imax,k1]=max(II);kmax=k1;
        %jet length (on-axis intensity above Imax/2)
        k2=kmax;while k1>1 && II(k1-1)>=Imax/2, k1=k1-1; end
        while k2<length(zz) && II(k2+1)>=Imax/2, k2=k2+1; end
        %transverse profile at the maximum
        if s==1
            P=PB_all(:,kmax)';xp=xs;
        else
            xp=linspace(-4*10^-6,4*10^-6,1601);
            r=sqrt(xp.^2+(zz(kmax)+R)^2);phi=atan2(xp,zz(kmax)+R);
            EP=exp(j*k*(zz(kmax)+R));
            for n=-nmax:nmax
                EP=EP+j^n*a(n+nmax+1)*besselh(n,1,k*r).*exp(j*n*phi);
            end
            P=abs(EP).^2;
        end
        i0=find(abs(xp)==min(abs(xp)),1);
        iL=find(P(1:i0)<Imax/2,1,'last');
        iR=i0-1+find(P(i0:end)<Imax/2,1,'first');
        xL=interp1(P(iL:iL+1),xp(iL:iL+1),Imax/2);
        xR=interp1(P(iR-1:iR),xp(iR-1:iR),Imax/2);
        FOM(c,s,:)=[zz(kmax)*10^6 Imax (xR-xL)/(Ld/nbg) (zz(k2)-zz(k1))*10^6];
        Prof{c,s}=[xp*10^6;P];
    end
    Axis{c,1}=[zB*10^6;IB];Axis{c,2}=[zM*10^6;IM];
    d_geo(c)=R*(2-m_rel(c))/(2*(m_rel(c)-1))*10^6;
end

%Results : BPM / Mie and relative error in %
fprintf('\n%-26s %6s | %15s | %15s | %15s | %15s\n','case','m','d [um]','eta=Imax/I0','FWHM/(lambda/n)','L1/2 [um]');
for c=1:length(np_list)
    fprintf('%-26s %6.3f |',names{c},m_rel(c));
    for q=1:4
        fprintf(' %6.2f / %6.2f |',FOM(c,1,q),FOM(c,2,q));
    end
    fprintf('\n');
end
Error_percent=100*(FOM(:,1,:)-FOM(:,2,:))./FOM(:,2,:);
Error_percent=squeeze(Error_percent)  %rows: cases, columns: d, eta, FWHM, L1/2

figure(1)
for c=1:length(np_list)
    subplot(2,2,c)
    plot(Axis{c,2}(1,:),Axis{c,2}(2,:),'k',Axis{c,1}(1,:),Axis{c,1}(2,:),'b--')
    grid on
    title(sprintf('%s, m = %.3f',names{c},m_rel(c)))
    xlabel('distance from the shadow surface [\mum]')
    ylabel('I(0,z) / I_0')
    legend('Mie (exact)','BPM')
end

figure(2)
for c=1:length(np_list)
    subplot(2,2,c)
    plot(Prof{c,2}(1,:),Prof{c,2}(2,:),'k',Prof{c,1}(1,:),Prof{c,1}(2,:),'b--')
    axis([-3 3 0 1.1*max([Prof{c,1}(2,:) Prof{c,2}(2,:)])])
    grid on
    title(sprintf('Profile at the jet maximum, m = %.3f',m_rel(c)))
    xlabel('x [\mum]')
    ylabel('I(x) / I_0')
    legend('Mie (exact)','BPM')
end

figure(3)
plot(m_rel,abs(Error_percent),'o-')
grid on
legend('d','\eta','FWHM','L_{1/2}','Location','northwest')
title('BPM error with respect to the exact Mie solution')
xlabel('relative index m = n_p / n_{bg}')
ylabel('|error| [%]')
