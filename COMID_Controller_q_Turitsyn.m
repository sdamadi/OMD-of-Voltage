function [q_t] = COMID_Controller_q_Turitsyn(g_q,q_g,PV_n,q_g_max,q_g_min,eta_q,c_n)

y_s = q_g - eta_q*g_q;

for i=1:PV_n
    
   if  q_g_max(i) + eta_q*c_n <= y_s(i)
    q_t(i) = q_g_max(i);
   
    elseif y_s(i) > eta_q*c_n & y_s(i)< q_g_max(i) + eta_q*c_n
    q_t(i) = y_s(i) - eta_q*c_n;
    
    elseif y_s(i) >= -eta_q*c_n & y_s(i)<= eta_q*c_n
    q_t(i) = 0;
    
    elseif y_s(i) >= q_g_min(i) - eta_q*c_n & y_s(i)< -eta_q*c_n
    q_t(i) = y_s(i) + eta_q*c_n;
    
    elseif y_s(i)<= q_g_min(i)- eta_q*c_n
    q_t(i) = q_g_min(i);  
          
end


end

