function [vm,delta,I_line] = PhaseRecovery_I(n,nbr,nf,nt,r,x,P_r,Q_r)
i=sqrt(-1);

V = zeros(n,1);

for k = 1:nbr
    
    V(1) = 1;
    
    V(nt(k)) = V(nf(k)) - ( r(k) + i*x(k) )*( P_r(k) - i*Q_r(k) )/V(nf(k))'; 
    
end

vm = abs(V);
delta = rad2deg(angle(V));

for k = 1:nbr     
 
    I_line(k,1) = ( P_r(k) - i*Q_r(k) )/V(nf(k))';
    
%     I_line_2(k,1) = ( V(nf(k)) - V(nt(k)) )/( r(k) + i*x(k) );
    
end



