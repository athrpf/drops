function qcf=qcfun(gl,ni,npt,dt,flag)
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
        izl=find((zd <20) | (zd >80));
        iyl=find((yd <10) | (yd >30));
        %iz2=find([1:npz]~=iz)
        qcfm(iyl,izl)=-3;
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);
    case 3
    % sinusfoermiger, stationaerer Waermestrom
        count=1;
        for i=0:dyd:yl
          qcfm(count,:)=sin((zd/zl)*4*pi);
          count=count+1;
        end
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);
    case 4
    % sinusfoermiger, instationaerer Waermestrom
    for t=1:npt
        count=1;
        for i=0:dyd:yl
            qcfm(count,:)=278.652e-3+200e-3*sin((zd/zl)*4*pi+(t/npt)*4*pi);
            %qcfm(count,:)=1e-3*sin((zd/zl)*2*pi+(t/npt)*2*pi);
            count=count+1;
        end
        qcfv=reshape(qcfm,npyz,1);
        qcf(:,t)=qcfv;
    end
    case 5
        qcfm=2*ones(npy,npz);
        izl=find((zd <20) | (zd >80));
        qcfm(:,izl)=0;
        iyl=find((yd<10) | (yd >30));
        qcfm(iyl,:)=0;
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);
    case 6
        qcfm=2*ones(npy,npz);
        qcfv=reshape(qcfm,npyz,1);
        qcf=repmat(qcfv,1,npt);
    case 7
    for t=1:npt
        count=1;
        for i=0:dyd:yl
            cnt=1;
            for j=0:dzd:zl
                qcfm(count,cnt)=i*(1-i)*j*(1-j)*(2+sin(t));
                %qcfm(count,cnt)=-2;
                cnt=cnt+1;
            end
            count=count+1;
        end
        qcfv=reshape(qcfm,npyz,1);
        qcf(:,t)=qcfv;
        %qcf=repmat(qcfv,1,npt);
    end
    case 8
    for t=1:npt
        count=1;
        for i=0:dyd:yl
            cnt=1;
            for j=0:dzd:zl
                qcfm(count,cnt)=2e3*i*(1-i)*j*(1-j)*sin(6*pi*t*dt);
                %qcfm(count,cnt)=-2;
%                 if (t<=50)
%                     qcfm(count,cnt)=100;
%                 else
%                     qcfm(count,cnt)=-100;
%                 end
                cnt=cnt+1;
            end
            count=count+1;
        end
        qcfv=reshape(qcfm,npyz,1);
        qcf(:,t)=qcfv;
        %qcf=repmat(qcfv,1,npt);
    end
end
