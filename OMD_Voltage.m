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
% ------------ % The variance of changing loads and PVs % ------------ %
var_1 = 0.1;
var_2 = 0.01;
% ------------ % The OMD's parameter % ------------ %
eta_q = 1;
eta_sig = 1;
c_til_0 = 6.6*1000;
c_n = 1/80;
c_til_n = 6.6*1000*c_n;
lambda = 10;
ep12 = 0;
ep13 = 0;
ep14 = 0;
ep15 = 0;
ep16 = 0.00;
ep23 = 0;
ep24 = 0;
ep25 = 0;
ep26 = 0.0;
ep34 = 0;
ep35 = 0;
ep36 = 0.0;
ep45 = 0;
ep46 = 0.1;
ep56 = 0.1;
    
C = [1,ep12,ep13,ep14,ep15,ep16;
     ep12,1,ep23,ep24,ep25,ep26;
     ep13,ep23,1,ep34,ep35,ep36;
     ep14,ep24,ep34,1,ep45,ep46;
     ep15,ep25,ep35,ep45,1,ep56;
     ep16,ep26,ep36,ep46,ep56,1];
Q = inv(C);
 
% ------------ % Number of periods % ------------ %
T = 121;
% ------------ % Number of realization of COMID % ------------ %
num_real = 1;
% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Preallocation of vectors to avoid wasting time % ------------ %
c0 = zeros(num_real,T);
c1 = zeros(num_real,T);
y2_s = zeros(n-1,1);
eps1 = zeros(n,1);
eps2 = zeros(n,1);
q_t = zeros(1,PV_n);
sigma_min_act = zeros(T,1);
f1 = zeros(T,1);
q_online = zeros(PV_n,1,T);
sigma_min_online = zeros(T,1);
p_g_rand = zeros(PV_n,1,T);
% sigma_min = zeros(T,1);
% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Initial values for running COMID % ------------ % 
q_g = zeros(PV_n,1);
q_g1 = zeros(PV_n,1);
q_g2 = zeros(PV_n,1);
sigma_min_0 = 0;
sigma_min = sigma_min_0;
sigma_min1 = sigma_min_0;
% ------------ ------------ ------------ ------------ ------------ %



% ------------ % Stochastic OMD for voltage stability % ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %
% ------------ % Stochastic OMD for voltage stability% ------------ %

% ------------ % Set the maximum and minimum q_g's that can be appled to the grid % ------------ %  

q_g_max = 0.45*p_g_nom;
q_g_min = -p_g_nom;




for o = 1:num_real
    
    
for k = 1:T
    
    k
    
    eta_q = 1/k;
    eta_sig = 1/k;

% ------------ % Set the power of PVs changing during the time % ------------ %    
    
p_g_rand(:,:,k) = p_g_nom + var_2*randn(PV_n,1);


% ------------ % Keep track of sigma_min evolution % ------------ %

sigma_min_online(k,1) = sigma_min;
sigma_min_online1(k,1) = sigma_min1;

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % The real loss of network (c0) without noise with Lambda % ------------ %

[c0(o,k),P_r,Q_r,l_r,vms_r,~,~,~] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand(:,:,k),q_g,children,injection_matrix,parent,PV_matrix,PV_n,lambda,zeros(n,1),zeros(n,1),0);
[c01(o,k),P_r1,Q_r1,l_r1,vms_r1,~,~,~] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand(:,:,k),q_g1,children,injection_matrix,parent,PV_matrix,PV_n,lambda,zeros(n,1),zeros(n,1),0);


% ------------ % The real loss of network (c0) without noise without Lambda % ------------ %

[c02(o,k),P_r2, Q_r2,l_r2,vms_r2] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand(:,:,k),q_g2,children,injection_matrix,parent,PV_matrix,PV_n,lambda,zeros(n,1),zeros(n,1),0);

% ------------ % Find the voltage pahses and current with Lambda% ------------ %

[~,delta,I_line] = PhaseRecovery_I(n,nbr,nf,nt,r,x,P_r,Q_r);
[~,delta1,I_line1] = PhaseRecovery_I(n,nbr,nf,nt,r,x,P_r1,Q_r1);

% ------------ % Find the voltage pahses and current without Lambda% ------------ %

[~,delta2,I_line2] = PhaseRecovery_I(n,nbr,nf,nt,r,x,P_r2,Q_r2);

% ------------ % Check the magnitude of current is the same as DistFlow % ------------ %

I_check = ( abs(I_line) ).^2 - l_r; 

% ------------ % Find the jacobian matrix and real minimum singular value of system with Lambda % ------------ %

cd 'C:\Users\Saeed\OneDrive\Matpower\matpower6.0'

filename='case47';
mpc=loadcase(filename);
mpc.bus(:,8) = sqrt(vms_r(2:end));
mpc.bus(:,9) = delta(2:end);
jac = full(makeJac(mpc));
sigma_min_act(k,1) = min(svd(jac));


filename='case47';
mpc=loadcase(filename);
mpc.bus(:,8) = sqrt(vms_r1(2:end));
mpc.bus(:,9) = delta1(2:end);
jac = full(makeJac(mpc));
sigma_min_act1(k,1) = min(svd(jac));

% ------------ % Find the jacobian matrix and real minimum singular value of system without Lambda % ------------ %

filename='case47';
mpc=loadcase(filename);
mpc.bus(:,8) = sqrt(vms_r2(2:end));
mpc.bus(:,9) = delta2(2:end);
jac = full(makeJac(mpc));
sigma_min_act2(k,1) = min(svd(jac));


cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\OMD-of-Voltage'

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Solve the problem with stochastic loads % ------------ %
% ------------ % Finding the dual variables y8_s corresponding to sigma_min from the dual problem with noise with Lambda % ------------ %

eps1(:,k) = randn(n,1);%zeros(n,1);
eps2(:,k) = randn(n,1);%zeros(n,1);

[c1(o,k),P_s,Q_s,l_s,vms_s,~,y2_s(:,k),~] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand(:,:,k),q_g,children,injection_matrix,parent,PV_matrix,PV_n,lambda,eps1(:,k),eps2(:,k),var_1);
[c11(o,k),P_s1,Q_s1,l_s1,vms_s1,~,y2_s1(:,k),~] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand(:,:,k),q_g1,children,injection_matrix,parent,PV_matrix,PV_n,lambda,eps1(:,k),eps2(:,k),var_1);


% ------------ % Finding the dual variables y8_s corresponding to sigma_min from the dual problem with noise without Lambda % ------------ %

[c12(o,k),P_s2,Q_s2,l_s2,vms_s2,~,y2_s2(:,k),~] = VoltageStability_DF_Solver_vms(nbr,n,r,x,p_c,q_c,p_g_rand(:,:,k),q_g2,children,injection_matrix,parent,PV_matrix,PV_n,lambda,eps1(:,k),eps2(:,k),var_1);


% ------------ % Cost regarding real loss c0 and q_g's achieved from OMD with Lambda % ------------ %

f1(k,1) = c_til_0*c0(o,k)+ c_til_n*sum(abs(q_g)) - lambda*sigma_min;   % c_til_0*c0(o,k) + c_til_n*sum(abs(q_g)) - lambda*sigma_min;  
f11(k,1) = c_til_0*c01(o,k)+ c_til_n*sum(abs(q_g1)) - lambda*sigma_min1;

% ------------ % Cost regarding real loss c0 and q_g's achieved from OMD without Lambda % ------------ %

f12(k,1) = c_til_0*c02(o,k) + c_til_n*sum(abs(q_g2));

% ------------ % The gradient of q_g (g_q) with Lambda % ------------ %
g_q = -y2_s([dual_indeces],k);
g_q_Test(:,1,k) = g_q;
g_q1 = -y2_s1([dual_indeces],k);

% ------------ % The gradient of q_g (g_q) without Lambda % ------------ %
g_q2 = -y2_s2([dual_indeces],k);

% ------------ % The gradient of sigma (g_sig) % ------------ %

if sigma_min_act(k,1) <= sigma_min  
    g_sig = 1;
elseif sigma_min_act(k,1) > sigma_min
    g_sig = -1;
else
end


if sigma_min_act1(k,1) <= sigma_min1  
    g_sig1 = 1;
elseif sigma_min_act1(k,1) > sigma_min1
    g_sig1 = -1;
else
end

% ------------ % Finding q_g for the next step with Lambda % ------------ %

C = inv(Q);

q_t = COMID_Controller_q(C,g_q,g_sig,q_g,PV_n,q_g_max,q_g_min,eta_q,c_n);

C = eye(size(Q));

q_t1 = COMID_Controller_q(C,g_q1,g_sig1,q_g1,PV_n,q_g_max,q_g_min,eta_q,c_n);



% ------------ % Finding q_g2 for the next step without Lambda % ------------ %

C = [eye(PV_n,PV_n),zeros(PV_n,1);zeros(1,PV_n),0];

q_t2 = COMID_Controller_q(C,g_q1,g_sig1,q_g2,PV_n,q_g_max,q_g_min,eta_q,c_n);


% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Finding sigma_min for the next step % ------------ %

C = inv(Q);
a = C*[g_q; g_sig];

if  0 < sigma_min - eta_q*a(6)
    sigma_min = sigma_min - eta_q*a(6);
   
    elseif 0 <= sigma_min - eta_q*a(6)
    sigma_min = 0;       
end

C = eye(size(Q));

a = C*[g_q1; g_sig1];

if  0 < sigma_min1 - eta_q*a(6)
    sigma_min1 = sigma_min1 - eta_q*a(6);
   
    elseif 0 <= sigma_min1 - eta_q*a(6)
    sigma_min1 = 0;       
end

%  sigma_min = sigma_min - eta_sig*g_sig;

% ------------ ------------ ------------ ------------ ------------ %


% ------------ % Keep track of all evolved q_g's with Lambda % ------------ %

q_g(:,1) = q_t;
q_online(:,1,k) = q_t;

q_g1(:,1) = q_t1;
q_online1(:,1,k) = q_t1;

% ------------ ------------ ------------ ------------ ------------ %

% ------------ % Keep track of all evolved q_g's without Lambda % ------------ %

q_g2(:,1) = q_t2;
q_online2(:,1,k) = q_t2;

% ------------ ------------ ------------ ------------ ------------ %

end
% ------------ % Check feasibility of the answer % ------------ % 
for i=1:nbr

    check_s(i,1) = norm ( [ 2*P_s(i) 2*Q_s(i) ( l_s(i) - vms_s(i) ) ] ) - ( l_s(i) + vms_s(i) ); 
    
end
% ------------ ------------ ------------ ------------ ------------ %

end


figure (1)
 
plot(0:T-1,f1,'--r',0:T-1,f11,'k')%,0:T-1,f12,':b'
xlabel('$t\,(min)$','Interpreter','latex')
xlim([0 T-1])
ylabel('$\tilde{c}_0 f_t(q^g)- \lambda\sigma_{COMID}$','Interpreter','latex')%ylabel('$\tilde{c}_0 f_t(q^g)+\tilde{c}_n\sum\limits_{n \in n_q }\left | q^g \right | - \lambda\sigma_{COMID}$','Interpreter','latex')
x_t = round(T/4,0);
y_t = f1(x_t,1);
% txt = ['$\sigma ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta_q,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f'),'$\,\,\,\,\sigma_{PV} ^ 2 =$',num2str(var_2,'%2.3f')];
txt = ['$\sigma ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta_q =$',num2str(eta_q,'%1.0f')];%txt = ['$\sigma ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta_q,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f')];
text(x_t,y_t,txt,'interpreter','latex')
legend({'$C=Q$','C=I'},'interpreter','latex')%,'$Withot \, \sigma$'

fig_1 = figure (1);
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\Figures_COMD-of-Voltage'
saveas(fig_1,sprintf('COMID_var_noise=%2.2f_eta=%1.0f_c_n=%2.5f_T=%d_var_PV=%2.3f.png',var_1,eta_q,c_n,T-1,var_2));
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\OMD-of-Voltage'

figure (2)
plot(0:T-1,sigma_min_act,'--r',0:T-1,sigma_min_act1,'k',0:T-1,sigma_min_act2,':b')
xlabel('$t\,(mins)$','Interpreter','latex')
xlim([0 T-1])
ylabel('$\sigma_{actual}$','Interpreter','latex')
x_t = round(T/4,0);
y_t = sigma_min_act(x_t,1);
% txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta_q,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f'),'$\,\,\,\,\sigma_{PV} ^ 2 =$',num2str(var_2,'%2.2f')];
txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta_q =$',num2str(eta_q,'%1.0f')];%txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta_q,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f')];

text(x_t,y_t,txt,'interpreter','latex')
legend({'$C=Q$','C=I','$Withot \, \sigma$'},'interpreter','latex')

fig_2 = figure (2);
cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\Figures_COMD-of-Voltage'
saveas(fig_2,sprintf('MinSingularValue_var_noise=%2.2f_eta=%1.0f_c_n=%2.5f_T=%d_var_PV=%2.3f.png',var_1,eta_q,c_n,T-1,var_2));

cd 'C:\Users\Saeed\OneDrive\UMBC\Dr. Kim\My papers\Matlab\First Paper\OMD-of-Voltage'


figure (3)
plot(0:T-1,sigma_min_online,'--r',0:T-1,sigma_min_online1,'k')
xlabel('$t\,(mins)$','Interpreter','latex')
xlim([0 T-1])
ylabel('$\sigma_{online}$','Interpreter','latex')
x_t = round(T/4,0);
y_t = sigma_min_act(x_t,1);
% txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta_q,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f'),'$\,\,\,\,\sigma_{PV} ^ 2 =$',num2str(var_2,'%2.2f')];
txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta_q =$',num2str(eta_q,'%1.0f')];%txt = ['$\sigma_{noise} ^ 2 =$',num2str(var_1,'%2.2f'),'$\,\,\,\,\eta =$',num2str(eta_q,'%1.0f'),'$\,\,\,\,c_n =$',num2str(c_n,'%2.5f')];%

text(x_t,y_t,txt,'interpreter','latex')
legend({'$C=Q$','C=I'},'interpreter','latex')


g1 = reshape(g_q_Test(1,1,1:T),1,T);
g2 = reshape(g_q_Test(2,1,1:T),1,T);
g3 = reshape(g_q_Test(3,1,1:T),1,T);
g4 = reshape(g_q_Test(4,1,1:T),1,T);
g5 = reshape(g_q_Test(5,1,1:T),1,T);

g_q_mean = [mean(g1);mean(g2);mean(g3);mean(g4);mean(g5)]
var_q_mean = [var(g1);var(g2);var(g3);var(g4);var(g5)]
nbins = 50;
 histogram(g1,nbins)
 histogram(g2,nbins)
 histogram(g3,nbins)
 histogram(g4,nbins)
 histogram(g5,nbins)
plot(2:n,sqrt(vms_r(2:n)),'-.or',2:n,sqrt(vms_r1(2:n)),'-.b')
plot(1:nbr,P_r(1:nbr),':b')
% powers available at each line must be associaated with its children, then
%  we have a chart based on buses gain.

[deltav,I] = max(abs(( sqrt(vms_r(3:n)) - vms_r(2)*ones(size(vms_r(3:n))) )/ vms_r(2) ));
[deltav1,I1] = max(abs(( sqrt(vms_r1(3:n)) - vms_r1(2)*ones(size(vms_r1(3:n))) )/ vms_r1(2) ));