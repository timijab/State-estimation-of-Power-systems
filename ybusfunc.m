%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN HELP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: ybus = ybusfunc(nbus) %%
%%%    %% 
%%% DESCRIPTION: The function ybusfun forms the admittance matrix %%
%%% of the IEEE bus system using resistance reactance %%
%%% and susceptance of the arrangement.  %%
%%% %%
%%% Arguments: nbus - Number of bus systems 14, 30 or 118 for IEEE %%
%%% 14 ,30 or 118 bus system. %%
%%%    %%
%%% Outputs: 1) ybus - The admittance matrix of the system.It %%
%%% is of the size n x n. Where n is number of bus. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END HELP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ybus = ybusfunc(nbus)  % Returns ybus
linedata = linedatas(nbus); % Calling "linedatas.m" for Line Data
fb = linedata(:,1);  % From bus 
tb = linedata(:,2);  % To bus 
r = linedata(:,3);  % Resistance, R
x = linedata(:,4);  % Reactance, X
b = linedata(:,5);  % half Ground Admittance, B/2
a = linedata(:,6);  % Tap setting. Its value is 1 for transmission line
z = r + 1i*x;  % Z matrix
y = 1./z;  % To get inverse of each element
b = 1i*b; 
busdata=busdatas(nbus);  % Calling "Busdatas.m" for shunt admittance
Bbus=1i*diag(busdata(:,3));
nbranch = length(fb);  % no. of branches
ybus = zeros(nbus,nbus);  % Initialize YBus
% Adding Off Diagonal Elements
for k=1:nbranch
ybus(fb(k),tb(k)) = ybus(fb(k),tb(k))-y(k)/a(k);
ybus(tb(k),fb(k)) = ybus(fb(k),tb(k));
end
% Adding Diagonal Elements
for m =1:nbus
ybus(m,m) = ybus(m,m)+Bbus(m,m);
for n =1:nbranch
if fb(n) == m
ybus(m,m) = ybus(m,m) + y(n)/(a(n)^2) + b(n);
elseif tb(n) == m
ybus(m,m) = ybus(m,m) + y(n) + b(n);
end
end
end