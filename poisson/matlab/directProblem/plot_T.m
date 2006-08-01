function plot_T(gl,ni,T,t,pp,flag)

yl=gl(2);
zl=gl(3);
niy=ni(2);
niz=ni(3);
npy=niy+1;
npz=niz+1;

% create mesh
dyd=yl/niy;
dzd=zl/niz;
yd=[0:dyd:yl];
zd=[0:dzd:zl];
[Yd,Zd]=meshgrid(yd,zd);

% create time-independent matrix
Tm=reshape(T(:,t),npy,npz);

if flag==1
    % set the limits of the x- and y-axis
    %xlim([0 yl]);
    %ylim([0 zl]);

    % output in the middle of the y-z-plane
    Tmp=Tm(niy/2+1,niz/2+1);

    if pp==1
        plot(t,Tmp,'linewidth',3);
    else
        plot(t,Tmp,'red--','linewidth',3);
    end
    %xlabel z
    %ylabel T
else
    if pp==1
        surfl(Yd,Zd,Tm');
        shading interp
    else
        mesh(Yd,Zd,Tm');
    end
    xlabel y
    ylabel z
    zlabel ('Temperatur T')
end


