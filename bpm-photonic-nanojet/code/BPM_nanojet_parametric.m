%BPM_nanojet_parametric.m
%Parametric study and improved design of the microcylinder nanojet (project part 2)
%Study 1 : relative index m = n_cylinder/n_water, D=10.3 um (BPM and exact Mie solution)
%Study 2 : cylinder size D/lambda, silica in water (BPM and exact Mie solution)
%Study 3 : improved design, silica core + low-index shell (design values n=1.35-1.42)
%Paramters common to all runs :
%Ld (light wavelength) =1.03, N=4096, L=160, dz=0.025 (micrometers)
%Run time : about 2-3 minutes (28 simulations). Results are printed run by run.
%Figures : 1 = study 1, 2 = study 2, 3 = choice of the design, 4 = intensity maps, 5 = on-axis intensity
%Each line of Runs is one simulation :
%[study  D(um)  n_core  n_shell  q=Rcore/R  nbg  w(um)  waist-centre(um)  behind(um)  Mie]
clear
%Gaussian Beam and grid
    N=4096; %N : sampling number
    L=160*10^-6; %Display area
    Ld=1.03;Ld=Ld*10^-6;
    k0=(2*pi)/Ld;  %wavenumber in vacuum, used for the index step
    dz=0.025;dz=dz*10^-6;
%Absorbing boundary (graded, starts at 80% of the half window)
    alpha=5*10^5;
%Silica and water at 1.03 um
    nSi=1.450;nW=1.326;
%List of simulations
    Runs=[];
    %Study 1 : relative index m (D=10.3 um in water)
    for m_rel=[1.05 1.10 1.15 1.20 1.30 1.45]  %stops at 1.45: BPM not reliable above (see validation)
        Runs=[Runs; 1 10.3 m_rel*nW m_rel*nW 1 nW 39.8 0 45 1];
    end
    %Study 2 : size D/lambda (silica in water)
    for Dl=[5 10 15 20 30]
        D=Dl*1.03;behind=max(45,1.2*(D/2)*(2-nSi/nW)/(2*(nSi/nW-1)));
        Runs=[Runs; 2 D nSi nSi 1 nW 39.8 0 behind 1];
    end
    %Study 3 : silica core + low-index shell (D=10.3 um), plus the homogeneous reference
    %shell indices 1.35-1.42 are design values (low-index polymer range)
    Runs=[Runs; 3 10.3 nSi nSi 1 nW 39.8 0 60 0];
    for ns=[1.35 1.37 1.40 1.42]
        for q=[0.4 0.5 0.6 0.7]
            Runs=[Runs; 3 10.3 nSi ns q nW 39.8 0 60 0];
        end
    end
%Useful jet : on-axis intensity above Imax/2 AND FWHM below Fmax*lambda/n
    Fmax=2;
%dx : step size
    dx=L/N;
    for m=1:N
     %Space axis
        x(m)=(m-1)*dx-L/2;
     %Frequency axis
        Kx(m)=(2*pi*(m-1))/(N*dx)-((2*pi*N)/(N*dx))/2;
    end
ic=N/2+1;
ix=find(abs(x)<=15*10^-6);xs=x(ix);i0=find(xs==0);
xa=0.8*L/2;
Abs=exp(-alpha*dz*(max(abs(x)-xa,0)/(L/2-xa)).^2);

for r=1:size(Runs,1)
    D=Runs(r,2)*10^-6;R=D/2;nc=Runs(r,3);ns=Runs(r,4);q=Runs(r,5);nbg=Runs(r,6);
    w0=Runs(r,7)*10^-6;off=Runs(r,8)*10^-6;behind=Runs(r,9)*10^-6;
    kbg=k0*nbg;
    zc=5*10^-6+R;Zw=zc+off;Z=zc+R+behind;
    %Gaussian Beam in space domain at z=0 (waist at z=Zw)
    Gau_w=exp(-x.^2./(w0^2));
    Gau_ini=ifft(fft(Gau_w).*fftshift(exp(j/(2*kbg)*(-Zw)*Kx.^2)));
    %Free space transfer function of propagation
    H=exp(j/(2*kbg)*dz*Kx.^2);
    %Iterative Loop
    Gau=Gau_ini;
    Nz=round(Z/dz);
    clear Zp I_axis P_pro
    for n1=1:Nz
        z=(n1-1)*dz;
        Zp(n1)=z+dz;
        %Gaussian Beam in Frequency domain
          FGau=fft(Gau);
        %Propagated Gaussian beam in Frequency domain
          FGau=FGau.*fftshift(H);
        %Propagated Gaussian beam in space domain
          Gau=ifft(FGau);
        %Refractive index map of the slice z+dz/2 (core of radius qR, shell up to R)
          rr=sqrt(x.^2+(z+dz/2-zc)^2);
          dn=(nc-nbg)*(rr<=q*R)+(ns-nbg)*(rr>q*R & rr<=R);
        %Step propagation through medium
          Gau=Gau.*exp(-j*k0*dn*dz).*Abs;
        %Stored intensity
          I_axis(n1)=abs(Gau(ic))^2;
          P_pro(:,n1)=abs(Gau(ix)').^2;
    end
    after=find(Zp>zc+R);
    zz=Zp(after)-(zc+R);II=I_axis(after);
    PB=P_pro(:,after);
    %Figures of merit of the BPM jet
    [Imax,kmax]=max(II);
    k1=kmax;k2=kmax;
    while k1>1 && II(k1-1)>=Imax/2, k1=k1-1; end
    while k2<length(II) && II(k2+1)>=Imax/2, k2=k2+1; end
    P=PB(:,kmax)';
    iL=find(P(1:i0)<P(i0)/2,1,'last');iR=i0-1+find(P(i0:end)<P(i0)/2,1,'first');
    if isempty(iL) || isempty(iR)  %no clean single peak (BPM failure at high contrast)
        xL=NaN;xR=NaN;
    else
        xL=interp1(P(iL:iL+1),xs(iL:iL+1),P(i0)/2);xR=interp1(P(iR-1:iR),xs(iR-1:iR),P(i0)/2);
    end
    lm=Ld/nbg;
    F_BPM(r,:)=[zz(kmax)*10^6 Imax (xR-xL)/lm (zz(k2)-zz(k1))/lm];
    %Useful jet length : longest segment with I > Imax/2 and FWHM(z) < Fmax*lambda/n
    good=zeros(1,length(zz));
    for t=1:4:length(zz)
        Pt=PB(:,t)';l=i0;h=i0;
        while l>1 && Pt(l-1)>=Pt(i0)/2, l=l-1; end
        while h<length(Pt) && Pt(h+1)>=Pt(i0)/2, h=h+1; end
        good(t:min(t+3,end))=(II(t)>=Imax/2)*((xs(h)-xs(l))/lm<=Fmax);
    end
    seg=diff([0 good 0]);b1=find(seg==1);b2=find(seg==-1)-1;
    if isempty(b1), L_use(r)=0; else, L_use(r)=max(zz(b2)-zz(b1))/lm; end
    Axis_r{r}=[zz*10^6;II];
    if Runs(r,1)==3
        Map_r{r}=P_pro(:,1:4:end);Zmap_r{r}=Zp(1:4:end)*10^6;Geo_r{r}=[zc R q]*10^6;
    end

    %Exact Mie solution (homogeneous cylinder, plane wave, E parallel to the axis)
    F_Mie(r,:)=NaN*[1 1 1 1];
    if Runs(r,10)==1
        m_rel=nc/nbg;k=kbg;X=k*R;Y=m_rel*X;nmax=round(X+4*X^(1/3)+15);
        clear a
        for n=-nmax:nmax
            Jx=besselj(n,X);Jy=besselj(n,Y);Hx=besselh(n,1,X);
            dJx=(besselj(n-1,X)-besselj(n+1,X))/2;dJy=(besselj(n-1,Y)-besselj(n+1,Y))/2;
            dHx=(besselh(n-1,1,X)-besselh(n+1,1,X))/2;
            a(n+nmax+1)=(m_rel*dJy*Jx-Jy*dJx)/(Jy*dHx-m_rel*dJy*Hx);
        end
        zM=linspace(R*1.0005,R+behind,3000);EM=exp(j*k*zM);
        for n=-nmax:nmax
            EM=EM+j^n*a(n+nmax+1)*besselh(n,1,k*zM);
        end
        IM=abs(EM).^2;[ImM,kM]=max(IM);
        k1=kM;k2=kM;
        while k1>1 && IM(k1-1)>=ImM/2, k1=k1-1; end
        while k2<length(IM) && IM(k2+1)>=ImM/2, k2=k2+1; end
        xp=linspace(-6*10^-6,6*10^-6,2401);rp=sqrt(xp.^2+zM(kM)^2);php=atan2(xp,zM(kM));
        EP=exp(j*k*zM(kM))*ones(size(xp));
        for n=-nmax:nmax
            EP=EP+j^n*a(n+nmax+1)*besselh(n,1,k*rp).*exp(j*n*php);
        end
        PP=abs(EP).^2;h=xp(PP>=ImM/2);
        F_Mie(r,:)=[(zM(kM)-R)*10^6 ImM (h(end)-h(1))/lm (zM(k2)-zM(k1))/lm];
    end
    fprintf('run %2d study %d | D=%5.1f nc=%5.3f ns=%5.3f q=%4.2f nbg=%5.3f w=%5.2f off=%5.0f | d=%6.2f eta=%5.2f FWHM=%4.2f L=%5.1f | Luse=%5.1f | Mie d=%6.2f eta=%5.2f FWHM=%4.2f L=%5.1f\n', ...
            r,Runs(r,1),Runs(r,2),nc,ns,q,nbg,Runs(r,7),Runs(r,8),F_BPM(r,:),L_use(r),F_Mie(r,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lab={'d [\mum]','\eta = I_{max}/I_0','FWHM / (\lambda/n)','L_{1/2} / (\lambda/n)'};

figure(1)  %Study 1 : relative index
s=find(Runs(:,1)==1);mm=Runs(s,3)./Runs(s,6);
for p=1:4
    subplot(2,2,p)
    plot(mm,F_Mie(s,p),'ko-',mm,F_BPM(s,p),'bs--')
    grid on
    xlabel('relative index m = n_p / n_{bg}')
    ylabel(lab{p})
    legend('Mie (exact)','BPM')
end
subplot(2,2,1);title('Study 1 : relative index (D = 10.3 \mum, water)')

figure(2)  %Study 2 : size
s=find(Runs(:,1)==2);DD=Runs(s,2)/1.03;
for p=1:4
    subplot(2,2,p)
    plot(DD,F_Mie(s,p),'ko-',DD,F_BPM(s,p),'bs--')
    grid on
    xlabel('D / \lambda_0')
    ylabel(lab{p})
    legend('Mie (exact)','BPM')
end
subplot(2,2,1);hold on
plot(DD,(DD*1.03/2)*(2-nSi/nW)/(2*(nSi/nW-1)),'r:')
hold off;legend('Mie (exact)','BPM','ray optics','Location','northwest')
title('Study 2 : size (silica in water)')

figure(3)  %Study 3 : useful jet length versus intensity
s=find(Runs(:,1)==3);base=s(1);
s1=find(Runs(:,1)==1);
plot(F_BPM(s1,2),L_use(s1),'b.-')
hold on
for ns=[1.35 1.37 1.40 1.42]
    c=s(2:end);c=c(Runs(c,4)==ns);
    plot(F_BPM(c,2),L_use(c),'o-')
end
plot(F_BPM(base,2),L_use(base),'k*','MarkerSize',12)
hold off
grid on
xlabel('\eta = I_{max} / I_0')
ylabel('useful jet length / (\lambda/n)')
legend('homogeneous cylinders (study 1)','silica core + shell 1.35','shell 1.37','shell 1.40','shell 1.42','silica cylinder (base)')
title(sprintf('Study 3 : jet length with I > I_{max}/2 and FWHM < %g \\lambda/n',Fmax))

%Improved design : longest useful jet with at least 90% of the base intensity (also printed: longest useful jet)
cand=s(2:end);
c1=cand(F_BPM(cand,2)>=0.9*F_BPM(base,2));[~,b]=max(L_use(c1));best=c1(b);
[~,b]=max(L_use(cand));longest=cand(b);
fprintf('\nBase       : silica cylinder D=10.3 um              | eta=%5.2f d=%5.2f um FWHM=%4.2f useful length=%5.1f lambda/n\n',F_BPM(base,2),F_BPM(base,1),F_BPM(base,3),L_use(base));
fprintf('Improved   : silica core q=%4.2f + shell n=%5.3f   | eta=%5.2f d=%5.2f um FWHM=%4.2f useful length=%5.1f lambda/n (x%4.2f)\n', ...
        Runs(best,5),Runs(best,4),F_BPM(best,2),F_BPM(best,1),F_BPM(best,3),L_use(best),L_use(best)/L_use(base));
fprintf('Longest    : silica core q=%4.2f + shell n=%5.3f   | eta=%5.2f d=%5.2f um FWHM=%4.2f useful length=%5.1f lambda/n (x%4.2f)\n', ...
        Runs(longest,5),Runs(longest,4),F_BPM(longest,2),F_BPM(longest,1),F_BPM(longest,3),L_use(longest),L_use(longest)/L_use(base));

figure(4)  %maps : base and improved design
th=linspace(0,2*pi,200);
tit={'Base : silica cylinder',sprintf('Improved : silica core (q = %.2f) + shell n = %.2f',Runs(best,5),Runs(best,4))};
rr=[base best];
for p=1:2
    subplot(2,1,p)
    imagesc(Zmap_r{rr(p)},xs*10^6,Map_r{rr(p)})
    axis xy;axis image
    ylabel(colorbar,'I / I_0')
    hold on
    g=Geo_r{rr(p)};
    plot(g(1)+g(2)*cos(th),g(2)*sin(th),'w')
    plot(g(1)+g(3)*g(2)*cos(th),g(3)*g(2)*sin(th),'w--')
    hold off
    axis([0 max(Zmap_r{rr(p)}) -10 10])
    title(tit{p})
    xlabel('z [\mum]');ylabel('x [\mum]')
end

figure(5)  %on-axis intensity : base and improved design
plot(Axis_r{base}(1,:),Axis_r{base}(2,:),'b',Axis_r{best}(1,:),Axis_r{best}(2,:),'r')
grid on
legend('base : silica cylinder','improved : core + shell')
xlabel('distance from the shadow surface [\mum]')
ylabel('I(0,z) / I_0')
title('On-axis intensity')