% Diese Funktion realisiert die graphische Ausgabe 
%---------------------------------------------------------

% Geometrieparameter
gl= [C_data(1),C_data(2),C_data(3)];

% Anzahl der Intervalle pro Koordinatenrichtung nach
% uniformer Verfeinerung
ni= [C_data(4),C_data(5),C_data(6)];

% nix= number of intervals in x-direction
% npx= number of points in x-direction
npx=ni(1)+1;
npy=ni(2)+1;
npz=ni(3)+1;
npyz= npy*npz;
dyd=gl(2)/ni(2);
dzd=gl(3)/ni(3);
dyz=1e-6*dyd*dzd;

% Zeitschrittweite
dt= C_data(7);
% Anzahl der Zeitschritte
ndt= C_data(8);

npt=ndt+1;
time=[0:dt:ndt*dt];

qcf=qcfun(gl,ni,npt,func);

%------------------- plot cost functional -----------------------

figure(1)
semilogy(J);xlabel('Iterations');ylabel('J');

return

%--------------- plot approximation over time -------------------

niter=size(J,2);
figure(2)
clear qcd fqc
for t=1:ndt+1
    qcd(:,:)=qc_iter(niter,:,:);
    plotq(gl,ni,qcd,t,0)
%     hold on
%     plotq(gl,ni,qcf,t,1)
%     hold off
    fqc(t)=getframe(gcf);
end

%-------------- plot approximation over iterations --------------

figure(3)
clear qcd fqc
time=60;
for i=1:niter
    qcd(:,:)=qc_iter(i,:,:);
    plotq(gl,ni,qcd,time,0)
    hold on
    plotq(gl,ni,qcf,time,1)
    hold off
    fqc(i)=getframe(gcf);
end

%--------------- plot residual over time -------------------

figure(4)
clear qcd fqc
for t=1:ndt+1
    qcd(:,:)=qc_iter(niter,:,:);
    res_q(:,:)=qcf(:,:)-qcd(:,:);
    plotq(gl,ni,res_q,t,1)
    fqc(t)=getframe(gcf);
end

