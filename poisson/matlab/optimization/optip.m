clear;

% gemeinsame Daten der direkten Probleme holen
[nu,xl,yl,zl,M,theta,cgtol,cgiter,Flag,BndRef,nis,ni,func]=getDirectData;

gl=[xl yl zl];
%

%
% nix= number of intervals in x-direction
% npx= number of points in x-direction
%
npx=ni(1)+1;
npy=ni(2)+1;
npz=ni(3)+1;
npyz= npy*npz;
dyd=yl/ni(2);
dzd=zl/ni(3);
dyz=1e-6*dyd*dzd
dt= 0.01;
ndt= 100;
npt=ndt+1;
time=[0:dt:ndt*dt];

% Create measurement data
Td=simdp;
qcf=qcfun(gl,ni,npt,func);
Timess=simip(qcf);
%figure(1)
%plot((Td+Timess)');

% Optimization
qc=zeros(npyz, ndt+1);
niter=0;
niter_max=10;
epsT=1e-6;
while 1    
    niter=niter+1 
    tic
    Tihat=simip(qc);
    comptime(niter,1)=toc;
    dThat=2*fliplr(Tihat-Timess);
    norm_err=norm(dThat)^2*dt*dyz;

    qc_iter(niter,:,:)=qc;
    That_iter(niter,:,:)=Td+Tihat;
    if (norm_err <epsT) | (niter >= niter_max)
        break
    else
        tic
        lambda=simipad(dThat);
        comptime(niter,2)=toc;
        Tad=fliplr(lambda);
        if niter==1;
            gamma=0;
            Pn=Tad;
        else        
            dTad=Tad-Tad_old;
            gamma=norm(Tad)/norm(Tad_old);  
            Pn=Tad+gamma*Pn;
        end        
        tic
        Tisens=simipsens(-Pn);
        comptime(niter,3)=toc;
        beta=-sum(diag((dThat*Tisens')))/norm(Tisens)^2;
        qc=qc-beta*Pn;
        Tad_old=Tad;
        J(niter)=norm_err
    end
    %figure(2);
    %plot((Td+Tihat)')
    %figure(3);
    %T_hat=Td+Tihat    
    %T_hat=Timess+Td    
    %plotq(gl,ni,qc,10,0)
    %hold on
    %plotq(gl,ni,qcf,10,1)
    %hold off

    %keyboard
    
    %save data qc_iter J ;
    
end

figure(1)
tfig=10;
plotq(gl,ni,qcf,tfig,1);
hold on
qcd(:,:)=qc_iter(niter,:,:);
plotq(gl,ni,qcd,tfig,0);
hold off

%Get the movie
%figure(4)
%clear qcd fqc
%tfig=10;
%for i=1:niter
%    qcd(:,:)=qc_iter(i,:,:);
%     plotq(gl,ni,qcd,tfig,0)
%     hold on
%    plotq(gl,ni,qcf,tfig,1)
%     hold off
%    fqc(i)=getframe(gcf);
%end
%movie(fqc)

% figure(5)
% clear qcd fTmess
% tfig=50;
% Tmess=Timess+Td;
% for i=1:niter
%    qcd(:,:)=That_iter(i,:,:);
%    plotq(gl,ni,qcd,tfig,0)
%    hold on
%    plotq(gl,ni,Tmess,tfig,1)
%    hold off
%    fTmess(i)=getframe(gcf);
% end
% movie(fTmess)

figure(6)
plot(J);xlabel('Iterations');ylabel('J')
%Plot solution
%figure(1);
%plot(time,Tihat');grid;
%plot(time,lambda');grid;
