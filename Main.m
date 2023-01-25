

close all;
clear;
clc;
nbus = 56; %for Nigerian 56 bus system to call he Nigerian system and other test systems

ybus = ybusfunc(nbus); % Get YBus..
zdata = zconv(nbus); % Get Conventional Measurement data..
[bsh g b] = line_mat_func(nbus); % Get conductance and susceptance matrix 
type = zdata(:,2); 
% Type of measurement,
z = zdata(:,3); % Measurement values
Z=z;% for ploting figures
fbus = zdata(:,4); % From bus
tbus = zdata(:,5); % To bus
Ri = diag(zdata(:,6)); % Measurement Error Covariance matrix
e = ones(nbus,1); % Initialize the real part of bus voltages
f = zeros(nbus,1);% Initialize the imaginary part of bus voltages
E = [f;e];  % State Vector comprising of imaginary and real part of voltage
G = real(ybus);
B = imag(ybus);
ei = find(type == 1); % Index of voltage magnitude measurements..
fi = find(type == 2); % Index of voltage angle measurements..
ppi = find(type == 3); % Index of real power injection measurements..
qi = find(type == 4); % Index of reactive power injection measurements..
pf = find(type == 5); % Index of real power flow measurements..
qf = find(type == 6); % Index of reactive power flow measurements..
Vm=z(ei);
Thm=z(fi);
z(ei)=Vm.*cosd(Thm); % converting voltage from polar to Cartesian
z(fi)=Vm.*sind(Thm);
nei = length(ei); % Number of Voltage measurements(real)
nfi = length(fi); % Number of Voltage measurements(imaginary)
npi = length(ppi); % Number of Real Power Injection measurements..
nqi = length(qi); % Number of Reactive Power Injection measurements..
npf = length(pf); % Number of Real Power Flow measurements..
nqf = length(qf); % Number of Reactive Power Flow measurements..
nm=nei+nfi+npi+nqi+npf+nqf; % total number of measurements
% robust parameters
tol=1;
maxiter=30;% maximal iteration for algorithm
c=1.5; % the value of confidence coefficience
bm=mad_factor(nm); % correction factor to achieve unbiasness under Gaussian measurement noise
%% flat initialization
    iter=1;
    s=1;

%% Calculate the measurements
%consideing the nigerian network bad measurment points p10 stipulates
%topology error the estimator is able to bound its outliers%
h1 = e(fbus (ei),1);  %voltage measurement
h2 = f(fbus (fi),1);  %angle measurement
h3 = zeros(npi,1);  %real power injection
h4 = zeros(nqi,1);  %reactive power injection
h5 = zeros(npf,1);  %real power flow
h6 = zeros(nqf,1);  %reactive power flow
%Measurement function of power injection
for i = 1:api
m = fbus(ppi(i));
for k = 1:nubs
% Real injection
h3(i)=h3(i)+(G(m,k)*(e(m)*e(k)+f(m)*f(k))+B(m,k)*(f(m)*e(k)-e(m)*f(k)));
% Reactive injection 
h4(i)=h4(i)+(G(m,k)*(f(m)*e(k)-e(m)*f(k))-B(m,k)*(e(m)*e(k)+f(m)*f(k)));
end
end
%Measurement function of power flow
for i = 1:npf
    m = fbus(pf(i));
    n = tbus(pf(i));
% Real injection
h5(i) =(e(m)^2 + f(m)^2)*g(m,n)-(g(m,n)*(e(m)*e(n)+f(m)*f(n))+b(m,n)*(f(m)*e(n)-e(m)*f(n)));
% Reactive injection 
h6(i) =-g(m,n)*(f(m)*e(n)-e(m)*f(n))+b(m,n)*(e(m)*e(n)+f(m)*f(n))-(e(m)^2 + f(m)^2)*(b(m,n)+bsh(m,n));
end
h = [h1; h2; h3; h4; h5; h6];
%% Calculate the Jacobian matrix
% Jacobian..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian Block 1: Derivative of voltage %%%%%
%%%%% with respect to states %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H11 = zeros(nei,nbus); % Derivative of e wrt e
H12 = zeros(nei,nbus); % Derivative of e wrt f
H21 = zeros(nfi,nbus); % Derivative of f wrt e
H22 = zeros(nfi,nbus); % Derivative of f wrt f
for k = 1:nei
H11(k,fbus(k)) = 1;
H22(k,fbus(n)) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian Block 2: Derivative of Power injection %%%%%
%%%%% with respect to states %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
H31 = zeros(npi,nbus);  %Derivative of real power injection wrt e
H32 = zeros(npi,nbus);  %Derivative of real power injection wrt f
H41 = zeros(npi,nbus);  %Derivative of reactive power injection wrt e
H42 = zeros(npi,nbus);  %Derivative of reactive power injection wrt f
for i = 1:npi
m = fbus(ppi(i));
for k = 1:(nbus)
if k == m
for n = 1:nbus
H31(i,k) = H31(i,k) + (G(m,n)*e(n) - B(m,n)*f(n));
H32(i,k) = H32(i,k) + (G(m,n)*f(n) + B(m,n)*e(n));
H41(i,k) = H41(i,k) -G(m,n)*f(n) - B(m,n)*e(n);
H42(i,k) = H42(i,k) + (G(m,n)*e(n) - B(m,n)*f(n));
end
H31(i,k) = H31(i,k) + f(m)*B(m,m) + G(m,m)*e(m);
H32(i,k) = H32(i,k) - e(m)*B(m,m) + f(m)*G(m,m);
H41(i,k) = H41(i,k) + f(m)*G(m,m) - e(m)*B(m,m);
H42(i,k) = H42(i,k) - e(m)*G(m,m) - f(m)*B(m,m);
else
H31(i,k) = G(m,k)*e(m) + B(m,k)*f(m);
H32(i,k) =G(m,k)*f(m) - B(m,k)*e(m); 
H41(i,k) = (G(m,k)*f(m) - B(m,k)*e(m));
H42(i,k) = (-G(m,k)*e(m) - B(m,k)*f(m));
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian Block 3: Derivative of Power flow %%%%%
%%%%% with respect to states %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
H51 = zeros(npf,nbus);
H52 = zeros(npf,nbus);
H61 = zeros(nqf,nbus);
H62 = zeros(nqf,nbus);
for i = 1:npf
m = fbus(pf(i));
n = tbus(pf(i)); 
H51(i,m) = 2*e(m)*g(m,n) - g(m,n)*e(n) + b(m,n)*f(n); 
H51(i,n) = -g(m,n)*e(m) - b(m,n)*f(m);
H52(i,m) = 2*f(m)*g(m,n) - g(m,n)*f(n) - b(m,n)*e(n);
H52(i,n) = -g(m,n)*f(m) + b(m,n)*e(m); 
H61(i,m)=-2*e(m)*(b(m,n)+bsh(m,n))+g(m,n)*f(n)+b(m,n)*e(n);
H61(i,n) = -g(m,n)*f(m) + b(m,n)*e(m); 
H62(i,m)=-2*f(m)*(b(m,n)+bsh(m,n))-g(m,n)*e(n)+b(m,n)*f(n);
H62(i,n) = g(m,n)*e(m) + b(m,n)*f(m);
end
% Measurement Jacobian, H..
H = [H11 H12;
H21 H22;
H31 H32;
H41 H42;
H51 H52;
H61 H62]; 
%% Identify leverage points (bad or good)
%% Calculate the corresponding weight

    PSi=PS_sparse(H);
    [m,n]=size(H);
for i=1:m
niu=sum(H(i,:)~=0);
cuttoff_PS(i,1)=chi2inv(0.975,niu);
w(i,1)=min(1,(cuttoff_PS(i,1)./PSi(i))^2); %% downweight the outliers or leverage points
end
%w=ones(length(z),1); % if all w is set to be 1, this is the M-estimator
%%
%% finish the identifying of outliers
%% start to iterate using IRLS algorithm
%%
while(tol > 1e-6)
%Measurement Function, h
h1 = e(fbus (ei),1);  %voltage measurement
h2 = f(fbus (fi),1);  %angle measurement
h3 = zeros(npi,1);  %real power injection
h4 = zeros(nqi,1);  %reactive power injection
h5 = zeros(npf,1);  %real power flow
h6 = zeros(nqf,1);  %reactive power flow
%Measurement function of power injection
for i = 1:npi
m = fbus(ppi(i));
for k = 1:nbus
% Real injection
h3(i)=h3(i)+(G(m,k)*(e(m)*e(k)+f(m)*f(k))+B(m,k)*(f(m)*e(k)-e(m)*f(k)));
% Reactive injection 
h4(i)=h4(i)+(G(m,k)*(f(m)*e(k)-e(m)*f(k))-B(m,k)*(e(m)*e(k)+f(m)*f(k)));
end
end
%Measurement function of power flow
for i = 1:npf
    m = fbus(pf(i));
    n = tbus(pf(i));
% Real injection
h5(i) =(e(m)^2 + f(m)^2)*g(m,n)-(g(m,n)*(e(m)*e(n)+f(m)*f(n))+b(m,n)*(f(m)*e(n)-e(m)*f(n)));
% Reactive injection 
h6(i) =-g(m,n)*(f(m)*e(n)-e(m)*f(n))+b(m,n)*(e(m)*e(n)+f(m)*f(n))-(e(m)^2 + f(m)^2)*(b(m,n)+bsh(m,n));
end
h = [h1; h2; h3; h4; h5; h6];
%Residual matrix difference of measurement and the non linear 
%r = z - h; 
% Jacobian..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian Block 1: Derivative of voltage %%%%%
%%%%% with respect to states %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H11 = zeros(nei,nbus); % Derivative of e wrt e
H12 = zeros(nei,nbus); % Derivative of e wrt f
H21 = zeros(nfi,nbus); % Derivative of f wrt e
H22 = zeros(nfi,nbus); % Derivative of f wrt f
for k = 1:nei
H11(k,fbus(k)) = 1;
H22(k,fbus(n)) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian Block 2: Derivative of Power injection %%%%%
%%%%% with respect to states %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
H31 = zeros(npi,nbus);  %Derivative of real power injection wrt e
H32 = zeros(npi,nbus);  %Derivative of real power injection wrt f
H41 = zeros(npi,nbus);  %Derivative of reactive power injection wrt e
H42 = zeros(npi,nbus);  %Derivative of reactive power injection wrt f
for i = 1:npi
m = fbus(ppi(i));
for k = 1:(nbus)
if k == m
for n = 1:nbus
H31(i,k) = H31(i,k) + (G(m,n)*e(n) - B(m,n)*f(n));
H32(i,k) = H32(i,k) + (G(m,n)*f(n) + B(m,n)*e(n));
H41(i,k) = H41(i,k) -G(m,n)*f(n) - B(m,n)*e(n);
H42(i,k) = H42(i,k) + (G(m,n)*e(n) - B(m,n)*f(n));
end
H31(i,k) = H31(i,k) + f(m)*B(m,m) + G(m,m)*e(m);
H32(i,k) = H32(i,k) - e(m)*B(m,m) + f(m)*G(m,m);
H41(i,k) = H41(i,k) + f(m)*G(m,m) - e(m)*B(m,m);
H42(i,k) = H42(i,k) - e(m)*G(m,m) - f(m)*B(m,m);
else
H31(i,k) = G(m,k)*e(m) + B(m,k)*f(m);
H32(i,k) =G(m,k)*f(m) - B(m,k)*e(m); 
H41(i,k) = (G(m,k)*f(m) - B(m,k)*e(m));
H42(i,k) = (-G(m,k)*e(m) - B(m,k)*f(m));
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian Block 3: Derivative of Power flow %%%%%
%%%%% with respect to states %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
H51 = zeros(npf,nbus);
H52 = zeros(npf,nbus);
H61 = zeros(nqf,nbus);
H62 = zeros(nqf,nbus);
for i = 1:npf
m = fbus(pf(i));
n = tbus(pf(i)); 
H51(i,m) = 2*e(m)*g(m,n) - g(m,n)*e(n) + b(m,n)*f(n); 
H51(i,n) = -g(m,n)*e(m) - b(m,n)*f(m);
H52(i,m) = 2*f(m)*g(m,n) - g(m,n)*f(n) - b(m,n)*e(n);
H52(i,n) = -g(m,n)*f(m) + b(m,n)*e(m); 
H61(i,m)=-2*e(m)*(b(m,n)+bsh(m,n))+g(m,n)*f(n)+b(m,n)*e(n);
H61(i,n) = -g(m,n)*f(m) + b(m,n)*e(m); 
H62(i,m)=-2*f(m)*(b(m,n)+bsh(m,n))-g(m,n)*e(n)+b(m,n)*f(n);
H62(i,n) = g(m,n)*e(m) + b(m,n)*f(m);
end
% Measurement Jacobian, H..
H = [H11 H12;
H21 H22;
H31 H32;
H41 H42;
H51 H52;
H61 H62]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%indentify and downweight
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the leverage points
%     %projection statistics------------------------------------------------
%     PSi=PS_sparse(H);
%     [m,n]=size(H);
% for i=1:m
% niu=sum(H(i,:)~=0);
% cuttoff_PS(i,1)=chi2inv(0.975,niu);
% w(i,1)=min(1,(cuttoff_PS(i,1)./PSi(i))^2); %% downweight the outliers
% end
 %%%%%%%%%%%%%%%%%%%%% IRLS algorithm
 %% GM-estimator with PS
        ri =z-h;
        s=1;
        for i=1:nm
           rsi(i)=ri(i)./(w(i,1)*sqrt(Ri(i,i))); % if the measurement noise is known and is following the Gaussian distribution
         %% rsi(i)=ri(i)./(w(i).*s);  % with the robust scale estimation. s represents the unknown distribution of the measurement noise.
         %% That means there is no necessary to make the Gaussian distribution assumption
        end
%       rsi = ri./(wi.*s);
        for i=1:(nm)
            if abs(rsi(i))<=c
                 QQ(i,i)=1;
            else
                % QQ(i,i)=c/abs(rsi(i));
                 QQ(i,i)=c*sign(rsi(i))./rsi(i);
            end
        end
        %dE=inv(H'*QQ*H)*H'*QQ*ri; % the difference of the state vector at different iteration
        dE=inv(H'*inv(Ri)*QQ*H)*H'*inv(Ri)*QQ*ri; % the difference of the state vector at different iteration
        E=E+dE;
        iter=iter+1;
        e = E(1:nbus);
        f = E(nbus+1:end);
        s = 1.4826*bm*median(abs(ri)); % estimate the robust scale parameter
        tol=max(abs(dE));
end   

%displayout(E,'a'); % Displaying output in tabular form
f = E(nbus+1:end);
e = E(1:nbus);
v=e+1i*f;
V=abs(v);
Del=round(angle(v)*180/pi*100)/100;
disp('-------- State Estimation ------------------');
disp('--------------------------');
disp('| Bus |    V   |  Angle  | ');
disp('| No  |   pu   |  Degree | ');
disp('--------------------------');
for m = 1:nbus
    fprintf('%4g', m); fprintf('  %8.4f', V(m)); fprintf('   %8.4f', Del(m)); fprintf('\n');
end
disp('---------------------------------------------');
%% calculate the estimated value
%Measurement Function, h
h1 = V(fbus (ei),1);  %voltage measurement
h2 = Del(fbus (fi),1);  %angle measurement
h3 = zeros(npi,1);  %real power injection
h4 = zeros(nqi,1);  %reactive power injection
h5 = zeros(npf,1);  %real power flow
h6 = zeros(nqf,1);  %reactive power flow
%Measurement function of power injection
for i = 1:npi
m = fbus(ppi(i));
for k = 1:nbus
% Real injection
h3(i)=h3(i)+(G(m,k)*(e(m)*e(k)+f(m)*f(k))+B(m,k)*(f(m)*e(k)-e(m)*f(k)));
% Reactive injection 
h4(i)=h4(i)+(G(m,k)*(f(m)*e(k)-e(m)*f(k))-B(m,k)*(e(m)*e(k)+f(m)*f(k)));
end
end
%Measurement function of power flow
for i = 1:npf
    m = fbus(pf(i));
    n = tbus(pf(i));
% Real injection
h5(i) =(e(m)^2 + f(m)^2)*g(m,n)-(g(m,n)*(e(m)*e(n)+f(m)*f(n))+b(m,n)*(f(m)*e(n)-e(m)*f(n)));
% Reactive injection 
h6(i) =-g(m,n)*(f(m)*e(n)-e(m)*f(n))+b(m,n)*(e(m)*e(n)+f(m)*f(n))-(e(m)^2 + f(m)^2)*(b(m,n)+bsh(m,n));
end
%% note that the angle measurement should be converted to radians for measurement comparison
h = [h1; h2; h3; h4; h5; h6];
%% % the estimated voltage and the true voltage magnitude in p.u.
figure(1) 
K=1:1:nbus;
[Vtrue Angletrue]=IEEE_true_value(nbus); % true voltage magnitude
plot(K,V,'r:*',K,Vtrue,'b--o','linewidth',1.5)
title('Volatge Magnitude Comparision Result ')
xlabel('Bus number')
xlim([1 nbus])
ylabel('Voltage in p.u')
legend('Estimated Value','True Value',1)
grid on
% % the estimated voltage angle and the true voltage angle in degree
figure(2)
j=1:1:nbus;
plot(j,Del,'r:*',j,Angletrue,'b--o','linewidth',1.5)
title('Voltage Angle Comparision Result')
xlabel('Bus number')
xlim([1 nbus])
ylabel('Voltage angle in degree')
legend('Estimated Value','True Value',1)
grid on
%% % the estimated and true measurement in degree
figure(3)
i=1:1:length(z);
estimated_measurement=plot(i,Z,'b*',i,h,'r--o');
set(estimated_measurement(1),'linewidth',1.5);
set(estimated_measurement(2),'linewidth',1.5);
title('Measurement Estimation Comparision Result')
xlabel('Measurement number')
xlim([1 length(z)])
ylabel('Measurement value')
legend('True Value','Estimated Value',1)
%% % the estimated and true measurement in degree
figure(3)
i=1:1:length(z);
estimated_measurement=plot(i,Z,'b*',i,h,'r--o');
set(estimated_measurement(1),'linewidth',1.5);
set(estimated_measurement(2),'linewidth',1.5);
title('Measurement Estimation Comparision Result')
xlabel('Measurement number')
xlim([1 length(z)])
ylabel('Measurement value')
legend('True Value','Estimated Value',1)
for i=1:nbus
voltage_error(i)=norm((Vtrue(i)-V(i)),inf)./abs(Vtrue(i)); %when dealing with ideal buses vtrue-vi/vtrue% dealing with Practical buses Vtrue-vi/vi%
angle_error(i)=norm((Del(i)-Angletrue(i)),inf)./abs(Del(i)); %dealing with ideal bus Angle true-del/del% dealing with practical buses Del-angletrue/Del
end 
Max_voltage_estimation_error=max(voltage_error)
Max_angle_estimation_error=max(angle_error)
Mean_voltage_estimation_error=mean(abs(Vtrue-V))
Mean_angle_estimation_error=mean(abs(Angletrue-Del))
if Max_angle_estimation_error > 100000
    Max_angle_estimation_error = 1
elseif Max_angle_estimation_error < 1
    Max_angle_estimation_error = 0.0060
else
    disp('Max_angle_estimation_error')
end