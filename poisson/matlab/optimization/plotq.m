%function plotq(gl,ni,qcv,t,pp)
function plotq(gl,ni,qcv,t,pp,flag)
%
yl=gl(2);
zl=gl(3);
niy=ni(2);
niz=ni(3);
npy=niy+1;
npz=niz+1;

% Create mesh
dyd=yl/niy;
dzd=zl/niz;
yd=[0:dyd:yl];
zd=[0:dzd:zl];
[Yd,Zd]=meshgrid(yd,zd);

% Create matrix qcm from vector qcv
qcm=reshape(qcv(:,t),npy,npz);

% flag=0 - plot von qc ueber der y-z-Ebene
% flag=1 - plot von qc fuer y=0
% flag=2 - plot von qc fuer z=0
%flag=0;

if flag==2
    qcm_2d(:,1)=qcm(:,1);
    if pp==1
        plot(yd,qcm_2d);
    else
        plot(yd,qcm_2d,'red--','linewidth',2);
    end
    xlabel y
    ylabel q_c
elseif flag==1
    % Ausgabe in der Mitte der Folie
    qcm_2d(1,:)=qcm(50,:);
    %qcm_2d(1,:)=qcm(1,:);
    % Ausgabe fuer y=0
    %qcm_2d(1,:)=qcm(8,:);
    if pp==1
        %ylim([-0.015 0.015]);
        plot(zd,qcm_2d,'linewidth',1.5);
    else
        plot(zd,qcm_2d,'red','linewidth',2);
        %axis([0 39 6.25e-3 6.6e-3])
    end
    %xlabel z
    %xlabel ('x[mm]','fontsize',12)
    %ylabel ('q_{foil}[kW/m^2]','fontsize',12)
else
    if pp==1
        %colormap(gray)
        surfl(Yd,Zd,qcm');

        shading interp
    else
        mesh(Yd,Zd,qcm');
        %axis([0 20 0 40 6.2e-3 6.7e-3])
        %axis([0 20 0 40 6e-3 7e-3])
    end
    %xlabel y
    %ylabel z
    %xlabel ('z[mm]','fontsize',12)
    %ylabel ('x[mm]','fontsize',12)
    %zlabel ('q_{foil}[kW/m^2]','fontsize',12)
    %title ('nopt=100','fontsize',16)
end
%[c,h]=contour(Yd,Zd,Tt,[-1:0.1:1]);
%chi=1:10:100;
%clabel(c(:,chi),h(chi));
%colorbar;