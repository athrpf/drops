function exflux = q(t)

% different benchmark heat fluxes

%if (0<t & t<=0.6)
%    exflux=t;
%elseif (0.6<t & t<=1.2)
%    exflux=1.2-t;
%else
%    exflux=0;
%end

%if (t<0.3)
%     exflux=0;
% elseif (t<0.4)
%     exflux=10*(t-0.3);
% elseif (t<0.8)
%     exflux=1;
% elseif (t<0.9)
%     exflux=1-10*(t-0.8);
% else
%     exflux=0;
% end

%exflux= 0;
%if (t>0.4 & t<=1)
%    exflux= 1;
%end

exflux= sin(2*pi*t/1.4);

return