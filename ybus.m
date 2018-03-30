%  Bus Admittance Matrix
%  Copyright (c)  1998 by H. Saadat.

function[Ybus] = ybus(zdata)
nl=zdata(:,1); nr=zdata(:,2); R=zdata(:,3); X=zdata(:,4);
nbr=length(zdata(:,1)); nbus = max(max(nl), max(nr));
Z = R + j*X;                            %branch impedance
y= ones(nbr,1)./Z;                     %branch admittance
Ybus=zeros(nbus,nbus);          % initialize Ybus to zero
for k = 1:nbr;   % formation of the off diagonal elements
    if nl(k) > 0 & nr(k) > 0
    Ybus(nl(k),nr(k)) = Ybus(nl(k),nr(k)) - y(k);
    Ybus(nr(k),nl(k)) = Ybus(nl(k),nr(k));
    end
end
for n = 1:nbus       % formation of the diagonal elements
    for k = 1:nbr
    if nl(k) == n | nr(k) == n
    Ybus(n,n) = Ybus(n,n) + y(k);
    else, end
    end
end