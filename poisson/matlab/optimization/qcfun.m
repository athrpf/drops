function qcf=qcfun(gl,ni,npt,flag)
yl=gl(2);
zl=gl(3);
niy=ni(2);
niz=ni(3);
npy=niy+1;
npz=niz+1;
npyz=npy*npz;
% Create mesh
dyd=yl/niy;
dzd=zl/niz;
yd=[0:dyd:yl];
zd=[0:dzd:zl];
qcfm=zeros(npy,npz);
switch flag
    case 1
        izl=find((zd <20) | (zd >80));
        %iz2=find([1:npz]~=iz)
        qcfm(:,izl)=-3;
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);        
    case 2
        izl=find((zd <20) | (zd >80))
        iyl=find((yd <10) | (yd >30));
        %iz2=find([1:npz]~=iz)
        qcfm(iyl,izl)=-3;
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);
    case 3
    % sinusfoermiger & stationaerer Waermestrom
        k=2;
        count=1;
        for i=0:dyd:yl
          qcfm(count,:)=sin((zd/zl)*2*k*pi);
          count=count+1;
        end  
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);
end
