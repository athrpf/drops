function plotq(gl,ni,qcv,t,pp)
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
% for i=1:npz
%     qcm(npy-i+1,:)=qcv(npy*(i-1)+1:i*npy,t)';
% end

if pp==1  
    %colormap(gray)  
    surfl(Yd,Zd,qcm');
    shading interp
else
    mesh(Yd,Zd,qcm');
end
    xlabel y
    ylabel z    
    zlabel q_c
%[c,h]=contour(Yd,Zd,Tt,[-1:0.1:1]);
%chi=1:10:100;
%clabel(c(:,chi),h(chi));
%colorbar;