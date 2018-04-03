clc
clear
close all
% 
%-----------------------------------------------%
% 
% 
%-----------------------------------------------%

MVA_max = [0 30	0	0	0	0	0	0	0	0	0	0.67	0.45	0	0.89	0	0.07	0	0.67	0	0	0.45	2.23	0	0	0.45	0.2	0	0.13	0.13	0.2	0.07	0.13	0.27	0.2	0	0.27	0	0.45	1.34	0.13	0.67	0.13	0	0.45	0.2	0.45	0]';

busdata = [
%n  code  p_c     q_c     p_g    q_g    shunt
1	3	0.0000	0.0000	0.0000	0.0000	0.0000
2	1	10.8000	8.1000	0.0000	0.0000	6
3	1	0.0000	0.0000	0.0000	0.0000	0.0000
4	1	0.0000	0.0000	0.0000	0.0000	1.2
5	1	0.0000	0.0000	0.0000	0.0000	0.0000
6	1	0.0000	0.0000	0.0000	0.0000	0.0000
7	1	0.0000	0.0000	0.0000	0.0000	0.0000
8	1	0.0000	0.0000	0.0000	0.0000	0.0000
9	1	0.0000	0.0000	0.0000	0.0000	0.0000
10	1	0.0000	0.0000	0.0000	0.0000	0.0000
11	1	0.0000	0.0000	0.0000	0.0000	0.0000
12	1	0.2412	0.1809	0.0000	0.0000	0.0000
13	1	0.1620	0.1215	0.0000	0.0000	0.0000
14	2	0.0000	0.0000	1.5 	0.0000	0.0000
15	1	0.3204	0.2403	0.0000	0.0000	0.0000
16	1	0.0000	0.0000	0.0000	0.0000	0.0000
17	1	0.0252	0.0189	0.0000	0.0000	0.0000
18	2	0.0000	0.0000	0.4 	0.0000	0.0000
19	1	0.2412	0.1809	0.0000	0.0000	0.0000
20	2	0.0000	0.0000	1.5 	0.0000	0.0000
21	1	0.0000	0.0000	0.0000	0.0000	0.0000
22	1	0.1620	0.1215	0.0000	0.0000	0.0000
23	1	0.8028	0.6021	0.0000	0.0000	0.0000
24	2	0.0000	0.0000	1   	0.0000	0.0000
25	2	0.0000	0.0000	2   	0.0000	0.0000
26	1	0.1620	0.1215	0.0000	0.0000	0.0000
27	1	0.0720	0.0540	0.0000	0.0000	0.0000
28	1	0.0000	0.0000	0.0000	0.0000	0.0000
29	1	0.0468	0.0351	0.0000	0.0000	0.0000
30	1	0.0468	0.0351	0.0000	0.0000	0.0000
31	1	0.0720	0.0540	0.0000	0.0000	0.0000
32	1	0.0252	0.0189	0.0000	0.0000	0.0000
33	1	0.0468	0.0351	0.0000	0.0000	0.0000
34	1	0.0972	0.0729	0.0000	0.0000	0.0000
35	1	0.0720	0.0540	0.0000	0.0000	0.0000
36	1	0.0000	0.0000	0.0000	0.0000	0.0000
37	1	0.0972	0.0729	0.0000	0.0000	0.0000
38	1	0.0000	0.0000	0.0000	0.0000	1.8
39	1	0.1620	0.1215	0.0000	0.0000	0.0000
40	1	0.4824	0.3618	0.0000	0.0000	0.0000
41	1	0.0468	0.0351	0.0000	0.0000	0.0000
42	1	0.2412	0.1809	0.0000	0.0000	0.0000
43	1	0.0468	0.0351	0.0000	0.0000	0.0000
44	1	0.0000	0.0000	0.0000	0.0000	0.0000
45	1	0.1620	0.1215	0.0000	0.0000	0.0000
46	1	0.0720	0.0576	0.0000	0.0000	0.0000
47	1	0.1620	0.1296	0.0000	0.0000	0.0000
48	1	0.0000	0.0000	0.0000	0.0000	1.8];

linedata= [
%nbr nf  nt         r         x
1	1	2	1.1e-5	1.1e-6
2	2	3	0.00170	0.00530
3	3	14	1.1e-5	1.1e-6
4	3	4	0.00020	0.00060
5	4	5	0.00030	0.00060
6	4	15	0.00060	0.00020
7	4	16	0.00140	0.00030
8	5	21	0.00220	0.00040
9	5	6	0.00070	0.00120
10	6	27	0.00040	0.00010
11	6	7	0.00010	0.00020
12	7	28	0.00110	0.00040
13	7	8	0.00020	0.00030
14	8	33	0.00050	0.00010
15	8	9	0.00010	0.00010
16	9	41	0.00030	0.00010
17	9	40	0.00160	0.00030
18	9	42	0.00070	0.00020
19	9	36	0.00050	0.00010
20	9	10	0.00020	0.00020
21	10	11	0.00010	0.00010
22	10	43	0.00100	0.00030
23	11	12	0.00070	0.00050
24	11	47	0.00150	0.00080
25	12	48	0.00020	0.00010
26	12	13	0.00050	0.00030
27	16	19	0.00030	0.00010
28	16	17	0.00070	0.00010
29	17	18	1.1e-5	1.1e-6
30	19	20	1.1e-5	1.1e-6
31	21	22	0.00080	0.00060
32	21	26	0.00140	0.00030
33	22	25	1.1e-5	1.1e-6
34	22	23	0.00130	0.00030
35	23	24	1.1e-5	1.1e-6
36	28	32	0.00030	0.00010
37	28	29	0.00070	0.00020
38	29	30	0.00070	0.00020
39	30	31	0.00040	0.00010
40	33	34	0.00030	0.00010
41	34	35	0.00020	0.00000
42	36	37	0.00050	0.00010
43	36	38	0.00050	0.00030
44	36	39	0.00070	0.00010
45	43	44	0.00040	0.00010
46	44	45	0.00040	0.00010
47	44	46	0.00040	0.00010];


% ------------ % Extract bus numbers % ------------ %

busnum = busdata(:,1);

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Extract bus types % ------------ %

buscode = busdata(:,2);

% #3 is slack bus
% #1 is PQ bus (Load)
% #2 is PV bus (Photovoltaic)

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Find number of branches % ------------ %

nbr = length(linedata(:,1));

% ------------ ------------ ------------ ------------ ------------ %

        
% ------------ % Find buses that sends energy (sending-end)  % ------------ %

nf  = linedata(:,2);

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Find buses that receive energy (receiving-end)  % ------------ %

nt  = linedata(:,3);

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Extract line impedances % ------------ %

r   = linedata(:,4);
x   = linedata(:,5);

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Find zdata from line data to calculate Y matrix % ------------ %

zdata1 = nf(2:end) - ones(nbr-1,1);
zdata2 = nt(2:end) - ones(nbr-1,1);
zdata3 = r(2:end);
zdata4 = x(2:end);
zdata = [zdata1 zdata2 zdata3 zdata4];
y_bim = sparse(ybus(zdata));

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Find the number of buses  % ------------ %

n = max(max(nf),max(nt));

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Set power factors for loads % ------------ %

pf = 0.8;

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Set ratio of loads out of maximum % ------------ %

load_ratio = 0.7;

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Set ratio of shunts out of maximum % ------------ %

Shunt_ratio = 0.6;

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Set ratio of PV generation out of maximum % ------------ %

PV_ratio = 0.6;

% ------------ ------------ ------------ ------------ ------------ %



% ------------ % Extract active loads % ------------ %

p_c = load_ratio*pf*MVA_max;

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Extract shunts % ------------ %

shunt = busdata(:,7);

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Forming reactive loads and considering shunts % ------------ %

% Since shunt capacitor is a generator of constant reactive power, considering 
% it as a negative value is reactive consumptions is the same

q_c = load_ratio*sqrt(1-pf^2)*MVA_max - Shunt_ratio*shunt;

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Set active power of PV generators % ------------ %

p_g = PV_ratio*busdata(:,5);

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Find the number of PV and PQ buses % ------------ %

PV_n = 0; % The initial number of PV buses
PQ_n = 0; % The initial number of PQ buses


for i = 1:n   
    % ------------ % Identifying which bus is a PV bus % ------------ %
    if buscode(i) == 2
    % ------------ % Count the number of PV buses % ------------ %
    PV_n = PV_n + 1;
    % ------------ % Find the indecies of constraints that q_g's contribute to them % ------------ %
    dual_indeces(PV_n) = find(nt==i);
    % ------------ % Identifying which bus is a PQ bus % ------------ %    
    elseif  buscode(i) == 1
    % ------------ % Count the number of PQ buses % ------------ %
    PQ_n = PQ_n + 1;
    
    end
end

% ------------ % Find the matrix of PV and PQ buses % ------------ %

% PV_matrix transforms q_g's (a vector with size PV_n) into a vector with size
% equal to the buses of the grid. Hence, the term (q_c - PV_matrix*q_g) is
% the negative injection vector.


PV_matrix = zeros(n,PV_n);

j1 = 0;


for i = 1:n
    
if buscode(i) == 2
    j1 = j1 + 1;
    PV_matrix(i,j1) = 1;
    
end
end

% ------------ % Find the matrix of Children and Parent buses % ------------ %

% Children matrix transforms vectors of active and reactive line flows (P and Q)
% into the set of summation of children of coresponding bus.

% Injection matrix assigns the bus quantity into line quantity.

% Parent matrix finds the parent of a specific bus


        children = zeros(nbr,nbr);
injection_matrix = zeros(nbr,n);
          parent = zeros(nbr,n);
       

for i=1:nbr
           for j=1:nbr
               if nf(j) == nt(i)
                   
% ------------ % Find the children buses in line #i % ------------ %

                  children(i,j) = 1;
                  
% ------------ % Find the children buses in line #i % ------------ %                  
                  
               end
               
% ------------ % Find the number of bus that injection must be considered in line #i % ------------ %

           injection_matrix(i,nt(i))= 1;
           
% ------------ % Find the number of bus that injection must be considered in line #i % ------------ %

           
% ------------ % Find the parent bus in line #i % ------------ %
           
           parent(i,nf(i)) = 1;
           
% ------------ % Find the parent bus in line #i % ------------ %

           end
end



% ------------ % Maximum active power produced by PV's % ------------ % 
p_g_max = [1.5 0.4 1.5 1 2]';
% ------------ % Nominal active power produced by PV's % ------------ %
p_g_nom = PV_ratio*p_g_max;
% ------------ % Maximum reactive power produced by inverters of PVs % ------------ % 
q_g_max = 0.45*p_g_max;
q_g_min = -q_g_max;
% ------------ % Initial reactive power for running OMD % ------------ % 
q_g = zeros(PV_n,1);
% ------------ % Number of periods % ------------ %
T = 5;
% ------------ % Number of realization of OMD % ------------ %
num_real = 1;
% ------------ % The variance of changing loads and PVs % ------------ %
var_1 = 0.0;
var_2 = 0.0;
% ------------ % The OMD's parameter % ------------ %
eta = 2;
c_til_0 = 6.6*1000;
c_n = 1/80;
c_til_n = 6.6*1000*c_n;
lambda = 1;
sigma_min_0 = 10;
jac0 = sigma_min_0*eye(2*(n-2) - PV_n,2*(n-2) - PV_n);


% ------------ % Preallocation of vectors to avoid wasting time % ------------ %
sigma_min_act = zeros(T,1);
% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Stochastic OMD for voltage stability % ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %


for o = 1:num_real
    
    jac = jac0;

for k = 1:T
    k

% ------------ % Set the power of PVs changing during the time % ------------ %    
    
p_g_rand = p_g_nom + var_2*randn(PV_n,1);

% ------------ % Set the maximum and minimum q_g's that can be appled to the grid % ------------ %  

q_g_max = 0.45*p_g_rand;
q_g_min = -p_g_rand;

% ------------ % The real loss of network (c0) without noise % ------------ %

[c0(o,k),P_r,Q_r,l_r,vms_r,~,~,~,~] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand,q_g,children,injection_matrix,parent,PV_matrix,PV_n,lambda,jac,zeros(n,1),zeros(n,1),0);

% ------------ % Find the voltage pahses and current % ------------ %

[~,delta,I_line] = PhaseRecovery_I(n,nbr,nf,nt,r,x,P_r,Q_r);

% ------------ % Check the magnitude of current is the same as DistFlow % ------------ %

I_check = ( abs(I_line) ).^2 - l_r; 

% ------------ % Find the jacobian matrix % ------------ %

cd 'C:\Users\Saeed\OneDrive\Matpower\matpower6.0'

filename='case47';
mpc=loadcase(filename);
mpc.bus(:,8) = sqrt(vms_r(2:end));
mpc.bus(:,9) = delta(2:end);
jac = full(makeJac(mpc));
sigma_min_act(k,1) = min(svd(jac));

cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\OMD-of-Voltage'

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Solve the problem with stochastic loads % ------------ %
% ------------ % Finding the dual variables y8_s corresponding to sigma_min from the dual problem with noise % ------------ %

eps1(:,k) = randn(n,1);%zeros(n,1);
eps2(:,k) = randn(n,1);%zeros(n,1);

[c1(o,k),P_s,Q_s,l_s,vms_s,~,y2_s(:,k),~,y8_s] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand,q_g,children,injection_matrix,parent,PV_matrix,PV_n,lambda,jac,eps1(:,k),eps2(:,k),var_1);



cvx_begin 

variables  sigma_min 
dual variables  y8
      
c5 = -lambda*sigma_min;
      minimize c5
      
         y8 : jac'*jac - sigma_min*eye(2*(n-2)-PV_n,2*(n-2)-PV_n) >= 0
         
%          sigma_min>= 0
 cvx_end



% ------------ % Cost regarding real loss c0 and q_g's achieved from OMD % ------------ %

f1(k,1) = c_til_0*c0(o,k) + c_til_n*sum(abs(q_g));   

% ------------ % The gradient % ------------ %
g_q = -y2_s([dual_indeces],k);

% g_sig = -y8_s(1,1);

% ------------ % Finding q_g for the next step % ------------ %
y_s_q = q_g - eta*g_q;

for i=1:PV_n

if  q_g_max(i) + eta*c_n < y_s_q(i)
    q_t(i) = q_g_max(i);
   
    elseif y_s_q(i) > eta*c_n && y_s_q(i)<= q_g_max(i) + eta*c_n
    q_t(i) = y_s_q(i) - eta*c_n;
    
    elseif y_s_q(i) >= -eta*c_n && y_s_q(i)<= eta*c_n
    q_t(i) = 0;
    
    elseif y_s_q(i) >= q_g_min(i) - eta*c_n && y_s_q(i)< -eta*c_n
    q_t(i) = y_s_q(i) + eta*c_n;
    
    elseif y_s_q(i)< q_g_min(i)- eta*c_n
    q_t(i) = q_g_min(i);
        
end

% ------------ % Finding sigma_min for the next step % ------------ %

% sigma_min_0 = sigma_min_0 - eta*g_sig;

% ------------ ------------ ------------ ------------ ------------ %

end

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Keeping track of all evolved q_g's % ------------ %

q_g(:,1) = q_t;

q_online(:,1,k) = q_t;

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Keeping track of sigma_min evolution % ------------ %

% sigma_min_online(k,1) = sigma_min_0;

% ------------ ------------ ------------ ------------ ------------ %


end
% ------------ % Check feasibility of the answer % ------------ % 
for i=1:nbr

    check_s(i,1) = norm ( [ 2*P_s(i) 2*Q_s(i) ( l_s(i) - vms_s(i) ) ] ) - ( l_s(i) + vms_s(i) ); 
    
end
% ------------ ------------ ------------ ------------ ------------ %

end


figure (1)
 
plot(0:k-1,f1)
xlabel('$t\,(min)$','Interpreter','latex')
xlim([0 T-1])
ylabel('$\tilde{c}_0 E[f_t(q^g)]+\tilde{c}_n\sum\limits_{n \in n_q }\left | q^g \right | $','Interpreter','latex')
x_t = round(T/4,0);
y_t = f1(x_t,1);
txt = ['$\sigma ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f'),'$\,\,\,\,\sigma_{PV} ^ 2 =$',num2str(var_2,'%2.3f')];
text(x_t,y_t,txt,'interpreter','latex')
legend({'COMD(Stochastic)'},'interpreter','latex')

fig_1 = figure (1);
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\Figures_COMD-of-Voltage'
saveas(fig_1,sprintf('COMD_var_noise=%2.2f_eta=%1.0f_c_n=%2.5f_T=%d_var_PV=%2.3f.png',var_1,eta,c_n,T-1,var_2));
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\OMD-of-Voltage'

figure (2)
plot(0:k-1,sigma_min_act)
xlabel('$t\,(mins)$','Interpreter','latex')
xlim([0 T-1])
ylabel('$\sigma_{min}$','Interpreter','latex')
x_t = round(T/4,0);
y_t = sigma_min_act(x_t,1);
txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f'),'$\,\,\,\,\sigma_{PV} ^ 2 =$',num2str(var_2,'%2.2f')];
text(x_t,y_t,txt,'interpreter','latex')
legend({'\,Minimum \,Singular \.Value \,of \,Jacobian'},'interpreter','latex')

fig_2 = figure (2);
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\Figures_COMD-of-Loss'
saveas(fig_2,sprintf('MinSingularValue_var_noise=%2.2f_eta=%1.0f_c_n=%2.5f_T=%d_var_PV=%2.3f.png',var_1,eta,c_n,T-1,var_2));
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\COMD'



% 
% 
% 
% 
% 
% 
% 
% 
% % ------------ % Deterministic OMD % ------------ %
% % ------------ % Deterministic OMD % ------------ %
% % ------------ % Deterministic OMD % ------------ %
% % ------------ % Deterministic OMD % ------------ %
% % ------------ % Deterministic OMD % ------------ %
% 
% 
% % ------------ % Initial reactive power for running OMD % ------------ % 
% q_g = zeros(PV_n,1);
% % ------------ % Number of periods % ------------ %
% 
% 
% for o = 1:num_real
% 
% for k = 1:T
%     k
% 
% [c_d(o,k),~,~,~,~,y2_d(:,k)] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g,q_g,children,injection_matrix,parent,PV_matrix,PV_n,lambda,zeros(n,1),zeros(n,1),0);
% f_d(k,1) = 6.6*1000*c_d(o,k) + 6.6*1000*c_n*sum(abs(q_g)); 
% 
% % ------------ % The gradient % ------------ %
% g_d = -y2_d([dual_indeces],k);
% 
% % ------------ % Finding q_g for the next step % ------------ %
% y_d = q_g - eta*g_d;
% 
% 
% for i=1:PV_n
% 
% if  q_g_max(i) + eta*c_n < y_d(i)
%     q_t(i) = q_g_max(i);
%    
%     elseif y_d(i) > eta*c_n & y_d(i)<= q_g_max(i) + eta*c_n
%     q_t(i) = y_d(i) - eta*c_n;
%     
%     elseif y_d(i) >= -eta*c_n & y_d(i)<= eta*c_n
%     q_t(i) = 0;
%     
%     elseif y_d(i) >= q_g_min(i) - eta*c_n & y_d(i)< -eta*c_n
%     q_t(i) = y_d(i) + eta*c_n;
%     
%     elseif y_d(i)< q_g_min(i)- eta*c_n
%     q_t(i) = q_g_min(i);
%         
% end
% end
% 
% q_g(:,1) = q_t;
% 
% q_det(:,1,k) = q_t;
% 
% 
% end
% end
% 
% 
% figure (1)
%  
% plot(0:k-1,f1,0:k-1,f_d)
% xlabel('$t\,(min)$','Interpreter','latex')
% xlim([0 T-1])
% ylabel('$\tilde{c}_0 E[f_t(q^g)]+\tilde{c}_n\sum\limits_{n \in n_q }\left | q^g \right | $','Interpreter','latex')
% x_t = round(T/4,0);
% y_t = f1(x_t,1);
% txt = ['$\sigma ^ 2 =$',num2str(var,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f')];
% text(x_t,y_t,txt,'interpreter','latex')
% legend({'OMD(Stochastic)','OMD(Deterministic)'},'interpreter','latex')
% 
% 
% 
% %%
% 
% 
% 
%  [c2,~, ~,~,~,~,~,q_g_2] = Optimal_DistFlowSolver_vms_q_g(nbr,n,PV_n,r,x,p_c,q_c,p_g,q_g,children,injection_matrix,parent,PV_matrix,q_g_max,zeros(n,1),zeros(n,1),0);
% 
%  [c3,~, ~,~,~,~,~,q_g_3] = Optimal_DistFlowSolver_vms(nbr,n,PV_n,r,x,p_c,q_c,p_g,q_g,children,injection_matrix,parent,PV_matrix,zeros(n,1),zeros(n,1),0);
% 
%  
%  f2 = 6.6*1000*c2+ 6.6*1000*c_n*sum(abs(q_g_2)); 
%  
%  f3 = 6.6*1000*c3 + 6.6*1000*c_n*sum(abs(q_g_3)); 
%  
% figure (2)
%  
% plot(0:k-1,f1,0:k-1,f_d,0:k-1,f2*ones(1,T),0:k-1,f3*ones(1,T))
% xlabel('$t\,(min)$','Interpreter','latex')
% xlim([0 T-1])
% ylabel('$\tilde{c}_0 E[f_t(q^g)]+\tilde{c}_n\sum\limits_{n \in n_q }\left | q^g \right | $','Interpreter','latex')
% x_t = round(T/4,0);
% y_t = f1(x_t,1);
% txt = ['$\sigma ^ 2 =$',num2str(var,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f')];
% text(x_t,y_t,txt,'interpreter','latex')
% legend({'OMD(Stochastic)','OMD(Deterministic)','Convex Problem','Unconstrained $q^g$'},'interpreter','latex')
% 
% figure (3)
% 
% plot(0:k-1,c0,0:k-1,c_d,0:k-1,c2*ones(1,T),0:k-1,c3*ones(1,T))
% xlabel('$t\,(min)$','Interpreter','latex')
% xlim([0 T-1])
% ylabel('$E[f_t(q^g)]$','Interpreter','latex')
% x_t = round(T/4,0);
% y_t = c0(x_t);
% txt = ['$\sigma ^ 2 =$',num2str(var,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f')];
% text(x_t,y_t,txt,'interpreter','latex')
% legend({'OMD(Stochastic)','OMD(Deterministic)','Convex Problem','Unconstrained $q^g$'},'interpreter','latex')
%  
% 



%%






