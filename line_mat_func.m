%%%  line_mat_func.m  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN HELP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: [bbus,g,b] = line_mat_func(nbus) %%
%%% %% 
%%% DESCRIPTION: The function line_mat_func returns shunt and series %%
%%% admittance of the IEEE bus system stance reactance %%
%%% and susceptance of the arrangement this is not same %%
%%% bus admittance matrix. %%
%%% %%
%%% Arguments: nbus - Number of bus systems 14, 30 or 118 for IEEE %%
%%% 14 ,30 or 118 bus system. %%
%%% %%
%%% Outputs: 1) bbus - The shunt admittance matrix. That is %%
%%% bbus(m,k) is shunt admittance of the line %%
%%% connecting bus m to bus k lumped at bus m %% 
%%% 2) g - Conductance of series admittance matrix. %%
%%% that is g(m,k) is series conductance of line %%
%%% connecting bus m and bus k %%
%%% 3) b - susceptance of series admittance matrix. %%
%%%    that is b(m,k) is series conductance of line %%
%%% connecting bus m and bus k %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END HELP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line Data for B-Bus (Shunt Admittance)Formation.
function [bbus,g,b] = line_mat_func(nbus) 
linedata = linedatas(nbus); % Calling "linedatas.m" for Line Data
fb = linedata(:,1);  % From bus 
tb = linedata(:,2);  % To bus 
r = linedata(:,3);  % Resistance, R
x = linedata(:,4);  % Reactance, X 
b_sh = linedata(:,5);  % half Ground Admittance, B/2
a = linedata(:,6);  % Tap setting. Its value is 1 for transmission line
z = r + 1i*x;  % Z matrix
y = 1./z;  % To get inverse of each element
nbranch = length(fb);  % no. of branches...
%bbus is shunt admittance matrix of the lines bbus (i,j) and bbus (j,i) are
%values of shunt admitance where i is sending end and j is receiving end.
%for line both are same but for transformer they are different.
%Y is series admittance matrix. g and b are series conductance and series
%susceptance respectively.
bbus = zeros(nbus);  %Initializing bbus
Y = zeros(nbus);  %Initializing Y
for k=1:nbranch
%Bbus for a Transmission line 
if a(k)==1
bbus(fb(k),tb(k)) = b_sh(k);
bbus(tb(k),fb(k)) = bbus(fb(k),tb(k));
%bbus for transformer with tap
else
bbus(fb(k),tb(k)) = imag(y(k))*(1-a(k))/a(k)^2;
bbus(tb(k),fb(k)) = imag(y(k))*(a(k)-1)/a(k);
end
% No condition is required for series admittance as tap ratio is 1 for 
% transmission line 
Y(fb(k),tb(k)) = y(k)/a(k);
Y(tb(k),fb(k)) = Y(fb(k),tb(k));
end
g=real(Y);
b=imag(Y);
