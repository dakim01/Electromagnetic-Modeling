
%BPM_nanojet_movie.m
%Movie of the propagation through the microcylinder
%Run BPM_nanojet_cylinder.m first (uses I_pro, I_axis, Zp, xs, zc, R, zs, kk, th,
%Imax, kmax, k1, k2, d_jet_um, d_geo_um, FWHM_um, xL, xR, L_half_um, Enhancement)
%Paramters of the movie :
%Nframe (number of frames) =160, FrameRate (frames per second) =20, Nhold (final frames) =40
%Run the simulation first if its results are not in the workspace
if ~exist('I_pro','var') || ~exist('kk','var')
    BPM_nanojet_cylinder
end
%Movie parameters
    %Nframe=input('number of frames = ');
    Nframe=160;
    %FrameRate=input('frames per second = ');
    FrameRate=20;
    Nhold=40;  %the last frame is held 2 s with the figures of merit
    step=floor(length(Zp)/Nframe);
    Imov=max(I_pro(:));  %fixed colour scale for all the frames
%Video file (MPEG-4 on Windows/Mac, AVI on Linux)
    if ispc || ismac
        v=VideoWriter('nanojet_cylinder.mp4','MPEG-4');
    else
        v=VideoWriter('nanojet_cylinder.avi','Motion JPEG AVI');
    end
    v.FrameRate=FrameRate;
    open(v);

%Field revealed slice by slice (zero where the beam has not arrived yet)
I_show=zeros(size(I_pro));
zs_um=zs*10^6;
frames=[step:step:length(Zp) length(Zp)*ones(1,Nhold)];

figure(6)  %figures 1-5 are used by BPM_nanojet_cylinder.m
set(gcf,'Color','w','Position',[100 100 1000 650])
for l=frames
    I_show(:,1:l)=I_pro(:,1:l);
    done=(l==length(Zp));
    %Intensity map |E(x,z)|^2 up to the current plane
      subplot(2,2,[1 2])
      imagesc(Zp,xs,I_show,[0 Imov])
      axis xy
      hold on
      plot((zc+R*cos(th))*10^6,R*sin(th)*10^6,'w')
      plot([Zp(l) Zp(l)],[-15 15],'w--')
      for q=1:4  %cross-section planes of task 6, drawn once the beam has passed them
          if kk(q)<=l, plot([Zp(kk(q)) Zp(kk(q))],[-15 15],'w:'), end
      end
      hold off
      axis image
      axis([min(Zp) max(Zp) -15 15])
      colorbar
      if done
          title(sprintf('Silica cylinder in water :  \\eta = %.1f,  d = %.1f \\mum,  FWHM = %.2f \\mum,  L_{1/2} = %.1f \\mum', ...
                Enhancement,d_jet_um,FWHM_um,L_half_um))
      else
          title(sprintf('Intensity |E(x,z)|^2 / I_0      z = %.1f \\mum',Zp(l)))
      end
      xlabel('z [\mum]')
      ylabel('x [\mum]')
    %Transverse profile at the current plane (FWHM drawn at the jet maximum)
      subplot(2,2,3)
      lp=l;if done, lp=kmax; end  %final frames : profile at the jet maximum
      plot(xs,I_pro(:,lp),'b')
      if done || abs(Zp(l)-Zp(kmax))<step*(Zp(2)-Zp(1))
          hold on
          plot([xL xR]*10^6,[Imax/2 Imax/2],'r-o')
          hold off
      end
      axis([-10 10 0 1.1*Imax])
      grid on
      title(sprintf('Transverse profile at z - z_s = %.1f \\mum',Zp(lp)-zs_um))
      xlabel('x [\mum]')
      ylabel('I(x,z) / I_0')
    %On-axis intensity up to the current plane
      subplot(2,2,4)
      plot(Zp(1:l),I_axis(1:l),'b')
      hold on
      plot([zs_um zs_um],[0 1.1*Imax],'k--')
      plot([zs_um zs_um]+d_geo_um,[0 1.1*Imax],'r:')
      if l>=kmax
          plot(Zp(kmax),Imax,'ro')
      end
      if l>=k2  %length L1/2 : on-axis intensity above Imax/2
          plot([Zp(k1) Zp(k2)],[Imax/2 Imax/2],'m-','LineWidth',2)
      end
      hold off
      axis([min(Zp) max(Zp) 0 1.1*Imax])
      grid on
      title('On-axis intensity')
      xlabel('z [\mum]')
      ylabel('I(0,z) / I_0')
    %Frame of the movie
      drawnow
      frame=getframe(gcf);
      writeVideo(v,frame);
end
close(v);
