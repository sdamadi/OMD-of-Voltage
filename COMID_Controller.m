function [q_t] = COMID_Controller(C,g_q,g_sig,q_g,PV_n,q_g_max,q_g_min,eta_q)
a = C*[g_q; g_sig];
for i=1:PV_n
    
        
if  q_g_max(i) <= q_g(i) - eta_q*a(i)
    q_t(i) = q_g_max(i);
   
    elseif q_g_max(i) > q_g(i) - eta_q*a(i) && q_g_min(i)< q_g(i) - eta_q*a(i)
    q_t(i) = q_g(i) - eta_q*a(i);
    
    elseif q_g_min(i) >= q_g(i) - eta_q*a(i)
    q_t(i) = q_g_min(i);    
          
end


end

