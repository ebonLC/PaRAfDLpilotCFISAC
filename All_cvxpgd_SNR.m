clc
clear all


%% Inital parameters
M=8; %number of access points
K=2; %number of terminals

D=1; %in kilometer

udp_ini = [0.2 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]; 
TAUu = 20; %Total pilot length

B=20e6; %Mhz
n_d = 10^((-164-30)/10);
noise_p = B*n_d; %= 10^((-203.975+10*log10(20*10^6)+9)/10);%
Hb = 15; % Base station height in m
Hm = 1.5; % Mobile height in m
f = 2400; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_f= 0.2; %downlink power: 200 mW
%noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
Pu = power_f/noise_p;%nomalized receive SNR
Pd=Pu;%pilot power

d0=0.01;%km
d1=0.05;%km

%Notice that in Hata-COST231 model for communication, the distance is
%measured in km while in radar it is in meter
AP = [0 0; 0 100; 100 0; 100 100; 200 0; 200 100; 300 0; 300 100]./1000; %AP coordinates

N=1000;
for nn = 1:N

nn

%Randomly locations of K terminals:
yk = 25:1:75; %[10:1:40 60:1:90];
UE_yk = yk(randperm(numel(yk),K)); %Randly generate y of UEs, continuous uniform distribution
xk = 50:1:250;
UE_xk = xk(randperm(numel(xk),K)); %Randly generate x of UEs
% sum(UE_x_k,2)/1000 sum(UE_y_k,2)/1000;
UE_ra(:,:,nn) = [UE_xk' UE_yk']./1000; 
% % UE(:,:,nn) = [0.1390 0.0450; 0.0520    0.0840];

%UE_rec(:,:,n) = [UE_x(n), UE_y; sum(UE_x_k,2)/1000, sum(UE_y_k,2)/1000]./1000;%UE_rec(:,:,n) = UE; %Record

%Test if the position generation is the same
for iii = 1:10
    for kk1 = 1:K
        for kk2 = 1:K
            Dist_test(kk1,kk2) = norm(UE_ra(kk1,:,nn) - UE_ra(kk2,:,nn));
            if Dist_test(kk1,kk2) == 0
                if ne(kk1,kk2)
                   UE_yk = yk(randperm(numel(yk),K)); %Randly generate y of UEs, continuous uniform distribution
                   UE_xk = xk(randperm(numel(xk),K));
                    %sum(UE_x_k,2)/1000 sum(UE_y_k,2)/1000;
                   UE_ra(:,:,nn) = [UE_xk' UE_yk']./1000;
                end
            end
        end
    end
end



UE = UE_ra(:,:,nn);


%Wrapped around (8 neighbor cells)
D = 1000/1000; %in 1000m...so in our case, virtual cells are not used

D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP(:,:,2) = AP(:,:,1) + D1; %right

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP(:,:,3)=AP(:,:,1)+D2; %up

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP(:,:,4)=AP(:,:,1)+D3; %left

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP(:,:,5)=AP(:,:,1)+D4; %down

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP(:,:,6)=AP(:,:,1)+D5; %right-down

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP(:,:,7)=AP(:,:,1)+D6; %left-up

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP(:,:,8)=AP(:,:,1)+D7; %right-up

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP(:,:,9)=AP(:,:,1)+D8; %left-down

%for n=1:N
%n

%UE = [UE_x(1), UE_y(1); UE_x(2), UE_y(2)]./1000; %UEs location
%UE = [150, UE_y(1); 150, UE_y(2)]./1000; %UEs location
%UE_rec(:,:,nn) = UE; %Record


sigma_shd=8; %in dB
D_cor=0.1;

% for z=1:length(udp_ini) %DL/UL ratio
% udp = udp_ini(z);
% T_u = TAUu*udp;%training length DL and sensing
% tau = TAUu - T_u;%training length UL

udp_ini = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]; %Pilot ratio
TAUu_total = 40;
for z=1:length(udp_ini)
TAUuu = TAUu_total*udp_ini(z); %The pilot use
T_u = round(TAUuu*0.5); %training length DL and sensing
tau = TAUuu - T_u; %training length UL


ii_oo = 1;
for iii = 1:ii_oo %times of random channel generation


%%%%%%Create the MxK correlated shadowing matrix %%%%%%%

    %%%%M correlated shadowing cofficients of M APs:
    Dist=zeros(M,M);%distance matrix
    Cor=zeros(M,M);%correlation matrix

    for m1=1:M
        for m2=1:M
            Dist(m1,m2) = min([norm(AP(m1,:,1)-AP(m2,:,1)), norm(AP(m1,:,1)-AP(m2,:,2)),norm(AP(m1,:,1)-AP(m2,:,3)),norm(AP(m1,:,1)-AP(m2,:,4)),norm(AP(m1,:,1)-AP(m2,:,5)),norm(AP(m1,:,1)-AP(m2,:,6)),norm(AP(m1,:,1)-AP(m2,:,7)),norm(AP(m1,:,1)-AP(m2,:,8)),norm(AP(m1,:,1)-AP(m2,:,9)) ]); %distance between AP m1 and AP m2
            Cor(m1,m2)=exp(-log(2)*Dist(m1,m2)/D_cor);
        end
    end
    A1 = chol(Cor,'lower');
    x1 = randn(M,1);
    sh_AP = A1*x1;
    for m=1:M
        sh_AP(m)=(1/sqrt(2))*sigma_shd*sh_AP(m)/norm(A1(m,:));
    end

    %%%%K correlated shadowing matrix of K terminal:
    Dist=zeros(K,K);%distance matrix
    Cor=zeros(K,K);%correlation matrix

    for k1=1:K
        for k2=1:K
            Dist(k1,k2)=norm(UE(k1,:) - UE(k2,:)); %distance between Terminal k1 and Terminal k2
            Cor(k1,k2)=exp(-log(2)*Dist(k1,k2)/D_cor);
        end
    end
    A2 = chol(Cor,'lower');
    x2 = randn(K,1);
    sh_Ter = A2*x2;
    for k=1:K
        sh_Ter(k)=(1/sqrt(2))*sigma_shd*sh_Ter(k)/norm(A2(k,:));
    end

%%% The shadowing matrix:
Z_shd=zeros(M,K);
for m=1:M
    for k=1:K
        Z_shd(m,k)= sh_AP(m)+ sh_Ter(k);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create an MxK large-scale coefficients beta_mk
BETAA = zeros(M,K);
dist=zeros(M,K);
for m=1:M  
    for k=1:K
        R(m,k) = norm(AP(m,:,1) - UE(k,:));
        [dist(m,k), index(m,k)] = min([ norm(AP(m,:,1) - UE(k,:)),norm(AP(m,:,2) - UE(k,:)),norm(AP(m,:,3) - UE(k,:)),norm(AP(m,:,4) - UE(k,:)),norm(AP(m,:,5) - UE(k,:)),norm(AP(m,:,6) - UE(k,:)),norm(AP(m,:,7) - UE(k,:)),norm(AP(m,:,8) - UE(k,:)),norm(AP(m,:,9) - UE(k,:)) ]); %distance between Terminal k and AP m
    if dist(m,k)<d0
         betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
    elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
         betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
    else
    betadB = -L - 35*log10(dist(m,k)) + Z_shd(m,k); %large-scale in dB
    end
    BETAA(m,k)=10^(betadB/10); 
    end

end

%% Collabrative Localization
h_los = ones(8,8).*0.85; %the vlaue of |RCS|^2

R_r = dist.*1000;
AP_r = AP.*1000;
UE_r = UE.*1000;

alpha = zeros(M,M,K); %Variation of signal strength
for i = 1:M
    for j = 1:M
        for k = 1:K
        alpha(i,j,k) = (3e8/(2.4*1e9))^2/(4*pi)^3 * 1/(R_r(i,k)^2 * R_r(j,k)^2); %Radar pathloss variation, carrier frequency = 1.9GHz
        end
    end
end

sigw = noise_p; %10^((-100-30)/10); % the noise power for target localization
c = 3e8; %Speed of light
Xi = 8*(pi^2)*(B^2)/(sigw*c^2); %Coefficient for simplification
for i = 1:M
    for j = 1:K
        xtr(i,j) = (AP_r(i,1,index(i,j))-UE_r(j,1))/R_r(i,j); %angle of x-axis
        ytr(i,j) = (AP_r(i,2,index(i,j))-UE_r(j,2))/R_r(i,j); %angle of y-axis
    end
end
ga = zeros(M,K); gb = zeros(M,K); gc = zeros(M,K);
for i = 1:M %Transmit AP
    for j = 1:M %Receive AP
        for k = 1:K %UE
            ga(i,j,k) = Xi*alpha(i,j,k)*(h_los(i,j))*(xtr(i,k)+xtr(j,k))^2;
            gb(i,j,k) = Xi*alpha(i,j,k)*(h_los(i,j))*(ytr(i,k)+ytr(j,k))^2;
            gc(i,j,k) = Xi*alpha(i,j,k)*(h_los(i,j))*(xtr(i,k)+xtr(j,k))*(ytr(i,k)+ytr(j,k));
        end
    end
end
Ga = sum(ga,2); Gb = sum(gb,2); Gc = sum(gc,2); %Each AP to all APs
GGa = reshape(Ga,M,K); GGb = reshape(Gb,M,K); GGc = reshape(Gc,M,K); %Reshape, column be index of transmit AP, row be index of UE
b_num = GGa + GGb;
A_den = zeros(M,M,K);
for i = 1:K
    A_den(:,:,i) = GGa(:,i)*GGb(:,i)' - GGc(:,i)*GGc(:,i)';
end

%% Pilot Asignment: (random choice) for UL
[U,S,V]=svd(randn(tau,tau));%U includes tau orthogonal sequences
Phii=zeros(tau,K);
% for k=1:K
%     Point=randi([1,tau]);
%     Phii(:,k)=U(:,Point);
% end

Point=randperm(tau,K); %If the number of UE is less than traning sequences, then no repeat
for k=1:K  
    Phii(:,k)=U(:,Point(k));
end

Phii_cf = Phii; % pilot set of cell-free systems
Phii_ULr = Phii_cf; %Record, used for conventional

% Create Gamma matrix (variances of the channel estimates)
Gammaa = zeros(M,K);
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau*Pu*BETAA(m,k)^2/(tau*Pu*mau(m,k) + 1);
    end
end

%When without DL, the UL estimation error is the error for the whole system
Var_mk003 = (BETAA - Gammaa) ./ BETAA; %normalized
Var_mk03 = sum(Var_mk003) ./ M; %system average to each UE
Var_mk3 = max(Var_mk03); %normalized worst average value

Etaa=zeros(M,1);
for m=1:M
    Etaa(m)=1/(sum(Gammaa(m,:))); %Water-filling
end
Pr_wop = ((noise_p)* T_u *Pd* (Etaa.*Gammaa));
P_u_full = sum(Pr_wop); 


%% Pilot Asignment: (random choice) for DL
[U0,S,V]=svd(randn(T_u,T_u));%U includes T_u orthogonal sequences
Phii12=zeros(T_u,K);

Point12=randperm(T_u,K);
for k=1:K  
    Phii12(:,k)=U0(:,Point12(k));
end

Phii_cf12 = Phii12; % pilot set of cell-free systems

thre = 1;

Phii_cf = Phii_cf12;

Phii_ULr = Phii_cf; %Record, used for conventional

%% Min-Max power optimizations for DL training + sensing - No pilot reuse

%tmax and tmin now should be the normalized value
tmin = min(min(Var_mk003))/100; 
tmax = Var_mk3; 

epsi = tmin/10;
PhiPhi=zeros(K,K);
for k=1:K
    for ii=1:K
        PhiPhi(k,ii)=norm(Phii_cf(:,k)'*Phii_cf(:,ii));
        PhiPhi2(k,ii)=PhiPhi(k,ii)^2;
    end
end

%% Full-Power 

%Full power localization MSE
Etaa=zeros(M,1);
for m=1:M
    Etaa(m)=1/(sum(Gammaa(m,:))); %Water-filling
end
Pr_wop = ((noise_p)* T_u *Pd* (Etaa.*Gammaa));
for q = 1:K
    b0(q) = b_num(:,q)'*Pr_wop(:,q);
    A0(q) = Pr_wop(:,q)'*A_den(:,:,q)*Pr_wop(:,q);
    Sig0(q) = b0(q)/A0(q);
end
L_mse_f = Sig0'; 
L_mse_full(iii) = ( L_mse_f(1) + L_mse_f(2) ) /2; 

%Full power channel estimation error
for k = 1:K %Order of the 'correct' UE
    for ii = 1:K %Order of all UE
        a_full(k,ii) = Etaa(:,1)'*(BETAA(:,k).*Gammaa(:,ii)) + PhiPhi2(k,ii).*( ( Etaa(:,1).^(1/2) )'*(BETAA(:,k)./BETAA(:,ii).*Gammaa(:,ii)) ).^2; %E{|channel|^2}
    end
end
for k = 1:K
        for kk = 1:K
            Skk_full(k,kk) = Etaa(:,1)'*(BETAA(:,k).*Gammaa(:,kk)); %Skk'
        end
    end
    for k = 1:K
        Var_minmax_full(k) = sum(Etaa(:,1).*BETAA(:,k).*Gammaa(:,k)) ./ (1+T_u*Pd* (Skk_full(k,:)*PhiPhi2(k,:)') ); %min-max error variance for each k
    end
    for k = 1:K
        nor_full(k) = Var_minmax_full(k)/a_full(k,k);
    end
[max_full index_full] = max(nor_full); %Finding the worst UE
Maxnor_full(iii) = max_full; 



 %% CVX within SCA+bisection
% 
Xopt = []; %initialize Xopt, for judging if the CVX is finally feasible
Xopt0 = [];

for m=1:M
    Xopt0(m)=1/(sum(Gammaa(m,:)));
end
Xopt = [Xopt0' Xopt0'];
Xoopt = Xopt;

Pr_0 = noise_p* T_u* Pd* (Xopt.*Gammaa); %Total power for DL
Pr00(iii) = sum(sum(Pr_0)); 

%cvx_solver sedumi
ITERAN = 20;
maxBis = 20; 

%Set a tnext, under the tnext, Inner finds the X which is close to tnext as
%possible, then Outer judges the tnext. If there's an X, lower than the
%bound, if not, raise the bound. 

for bis = 1:maxBis %Max bisection iterate times
    tnext = (tmax + tmin) / 2;
    Xopt_i = Xopt;
    
    for iter = 1:ITERAN
        %The initial F,W,Z and H
        for k = 1:K
            F_i(k) = Xopt_i(:,k)'*(BETAA(:,k).*Gammaa(:,k));
            Zsqr_i(k) = (Xopt_i(:,k).^(1/2))'*Gammaa(:,k); 
            Z_i(k) = ( Zsqr_i(k) )^2 + F_i(k);   
        end
        Wk_0 = zeros(K,1);
        for k = 1:K
            Wk_c(k) = 1;
            for ii = 1:K
                Skk(k,ii) = sum( Xopt_i(:,ii).* BETAA(:,k).* Gammaa(:,ii) ); %The Skk and Skk', variance of effective downlink channel
            end
            Wk_0(k) = Wk_c(k) + T_u*Pd* ( Skk(k,:) * PhiPhi2(k,:)' );
        end
        for k = 1:K
            H_i(k) = tnext* Wk_0(k)* Z_i(k) - F_i(k);
        end

        %Derivative of H_i
        Hde_i = zeros(M,K,K); %W(X)*Z'(X) + W'(X)*Z(X)
        dZ = zeros(M,K); 
        dW = zeros(M,K); 
        for k = 1:K
            dZ(:,k) = Zsqr_i(k)* (Gammaa(:,k)./( Xopt_i(:,k).^(1/2) ) )+ (Gammaa(:,k).*BETAA(:,k));
            for ii = 1:K
                dW(:,ii) = T_u*Pd* ( BETAA(:,k).*Gammaa(:,ii) )* PhiPhi2(k,ii);
            end
            for m = 1:M
                for j = 1:K %Since W depends on k'(j), we need to derived H by X_{m,j}(X_{m,k'}), not the X_{m,k}
                    Hde_i(m,j,k) = tnext* ( Z_i(k)* dW(m,j) + Wk_0(k)* dZ(m,j) ) - (j==k)*( Gammaa(m,k)* BETAA(m,k) );
                end
            end
        end

        %CVX, Taylor expansion of H should not less than zero
        cvx_begin quiet
            variable X(M,K) nonnegative %Constraint (28b)
            subject to
                for k = 1:K
                    H_tid = H_i(k)+ sum( sum( Hde_i(:,:,k).*( X - Xopt_i ) ) ); 
                    H_tid >= 0;
                end
                for m=1:M
                    X(m,:)* Gammaa(m,:)' <= 1; %Power constraint for each AP
                end
                Pr = (noise_p)* T_u*Pd* X.*Gammaa ; %normalized Pilot Power for corresponding UE
                for p = 1:K
                    b_num(:,p) - thre*A_den(:,:,p)*Pr(:,p) <= zeros(M,1);
                end
                minimize(0)
        cvx_end

        %Avoid inner is infeasible, which gives Xopt_i NaN
        if strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed')
            tmin = tnext;
            break; %Break out inner SCA loop immediately
        end
        
        Xopt_i = X; 
    
    end

    % bisection
    for k = 1:K
        Fs(k) = X(:,k)'*(BETAA(:,k).*Gammaa(:,k));
        Zsqr_s(k) = (X(:,k).^(1/2))'*Gammaa(:,k); 
        Zs(k) = ( Zsqr_s(k) )^2 + Fs(k);
    end
    Wk_s = zeros(K,1);
    for k = 1:K
        Wk_cs(k) = 1;
        for ii = 1:K
            Skk_s(k,ii) = sum( X(:,ii).* BETAA(:,k).* Gammaa(:,ii) ); %The Skk and Skk', variance of effective downlink channel
        end
        Wk_s(k) = Wk_cs(k) + T_u*Pd* ( Skk_s(k,:) * PhiPhi2(k,:)' );
    end
    A_s = Fs./ (Wk_s'.* Zs); 
    if max(A_s) <= tnext % feasible
    fprintf(1,'Problem is feasible ',tnext);
        tmax = tnext; 
        Xopt = X; 
    else % not feasible
    fprintf(1,'Problem not feasible ',tnext);
        tmin = tnext;   
    end

    if tmax - tmin < epsi
        break;
    end
end

if  strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed') %if the problem is infeasible, use full-power
    for k = 1:K %Order of the 'correct' UE
        for ii = 1:K %Order of all UE
            a_kr(k,ii) = Xoopt(:,ii)'*(BETAA(:,k).*Gammaa(:,ii)) + PhiPhi2(k,ii).*( ( Xoopt(:,ii).^(1/2) )'*(BETAA(:,k)./BETAA(:,ii).*Gammaa(:,ii)) ).^2; %E{|channel|^2}
        end
    end
    Pr = ((noise_p)* T_u*Pd* (Xoopt.*Gammaa));
    for k = 1:K
        for kk = 1:K
            Skkk(k,kk) = Xoopt(:,kk)'*(BETAA(:,k).*Gammaa(:,kk)); %Skk'
        end
    end
    for k = 1:K
        Var_minmax(k) = sum(Xoopt(:,k).*BETAA(:,k).*Gammaa(:,k)) ./ (1+T_u*Pd* (Skkk(k,:)*PhiPhi2(k,:)') ); %min-max variance for each k
    end
    for k = 1:K
        nor_Var(k) = Var_minmax(k)/a_kr(k,k);
    end
    Maxnor_Var(iii) = max(nor_Var); 
else %If problem has solution, just tnext is the answer, thenjust calculate Pr
    Pr = ((noise_p)* T_u*Pd* (Xopt.*Gammaa));
    for k = 1:K
        for kk = 1:K
            Skkk(k,kk) = Xopt(:,kk)'*(BETAA(:,k).*Gammaa(:,kk)); %Skk'
        end
    end
    for k = 1:K
        Var_minmax(k) = sum(Xopt(:,k).*BETAA(:,k).*Gammaa(:,k)) ./ (1+T_u*Pd* (Skkk(k,:)*PhiPhi2(k,:)') ); %min-max variance for each k
    end
    for k = 1:K %Order of the 'correct' UE
        for ii = 1:K %Order of all UE
            a_kr(k,ii) = Xopt(:,ii)'*(BETAA(:,k).*Gammaa(:,ii)) + PhiPhi2(k,ii).*( ( Xopt(:,ii).^(1/2) )'*(BETAA(:,k)./BETAA(:,ii).*Gammaa(:,ii)) ).^2; %E{|channel|^2}
        end
    end
    for k = 1:K
        nor_Var(k) = Var_minmax(k)/a_kr(k,k);
    end
    Maxnor_Var(iii) = tnext; 
end


%% Min-Max power allocations for No DL + sensing
% %Maintance the same radar SNR
% %Radar SINR
% for i = 1:M %Transmit AP
%     for j = 1:M %Receive AP
%         for k = 1:K %User
%             SNR(i,j,k) = (alpha(i,j,k)* h_los(i,j)* Pr(i,k) ) / (noise_p);
%         end
%     end
% end
% SNR_aver0 = sum(sum(SNR));%average radar SNR for whole system to each UE
% SNR_aver = reshape(SNR_aver0,[2,1]);
% 
% [Xopt,L_mse00,KK,Gammaaa,Pr_rr] = Without_DL_r(SNR_aver,M,K,Phii_ULr,tau,Pu,Pd,BETAA,Gammaa,b_num,A_den,noise_p,alpha,h_los,T_u)
% %max the communication performance while radar SNR ==
% %Assume power for data transmission is the same as DL, data power \-/ pilot power/symbol
% 
% L_mse000(iii) = ( L_mse00(1)+L_mse00(2) )/2;
% 
% Pr_rrr(iii) = sum(sum( Pr_rr ));

 %% PGD within feasible projection
% 
% %Avoid intereference
% gamma = Gammaa; 
% beta = BETAA; 
% C = PhiPhi2; 
% Desc   = 1e-15; %Gradient step size
% maxite = 1000;
% Fcc = noise_p * T_u * Pd;
% 
% X = ones(M,K) * (1/mean(gamma(:))); %Reset X
% X = max(X,0);
% for m = 1:M
%     s = gamma(m,:)*X(m,:)';
%     if s > 1
%         X(m,:) = X(m,:) * (1/s); %(28c)
%     end
% end
% 
% for it = 1:maxite
%     %F, W, Z as CVX
%     F = sum(beta.*gamma.*X)';  
%     W = zeros(K,1);
%     for k = 1:K
%         cosnt = 0;
%         for ii = 1:K
%             cosnt = cosnt + C(k,ii)*sum(X(:,ii).*beta(:,k).*gamma(:,ii));
%         end
%         W(k) = 1 + T_u*Pd * cosnt;
%     end
%     zsqr = sum(gamma.*sqrt(X),1)';
%     Z = zsqr.^2 + F;
%     A = F ./ (W .* Z);               
% 
%     %Worst A, the same as derivation, but notice gradiant A
%     [~, kmax] = max(A);
%     dF    = beta(:,kmax).*gamma(:,kmax);
%     dzsqr    = gamma(:,kmax)./(2*sqrt(X(:,kmax)   +  1e-12 ));
%     dZ    = 2*zsqr(kmax)*dzsqr + dF;
%     gradA = zeros(M,K);
%     gradA(:,kmax) = ( dF*Z(kmax) - F(kmax)*dZ ) / ( W(kmax)*Z(kmax)^2 ); %Quotient rule, dy/dx = ..../v^2
% 
%     %Descent step
%     Xn = X - Desc*gradA; 
%     Xn = real(Xn);
%     Xn = max(Xn, 0);
% 
%     %Constraint/ feasible projection
% for Rroww = 1: 15 %Xn on each projection, try to get the Xn close to three at the same time    
%     Xn = max(Xn, 0); %(28b)
% 
%     for p = 1:K
%         y = Xn(:,p).*gamma(:,p);
%         for i = 1:M
%             a     = thre * Fcc * A_den(i,:,p)';
%             viol  = b_num(i,p) - a'*y;
%             if viol > 0
%                 y = y + (viol/(a'*a)) * a; %(28d)
%             end
%         end
%         Xn(:,p) = y ./ gamma(:,p);
%     end
% 
%     for m = 1:M
%         u = Xn(m,:)';  g = gamma(m,:)';
%         s = g'*u;
%         if s > 1
%             u = u - ((s-1)/(g'*g))*g; %(28c), halfspace
%             
%             u_s = min(u);
%             if u_s < 0
%                 for m = 1:M
%                     s_s = gamma(m,:)*Xn(m,:)';
%                     if s_s > 1
%                         Xn(m,:) = Xn(m,:) * (1/s_s); %Avoid minus X and over power
%                         u = Xn(m,:)'; 
%                     end
%                 end
%             end
%         end
%         Xn(m,:) = u';
%     end
% end
% %If the iteration still cannot avoid negative or complex X...
% Xn = real(Xn);
% Xn = max(Xn,0); 
% 
%     if norm(Xn - X,'fro')/norm(X,'fro') < epsi/100
%         X = Xn;
%         break; %No need to reach the full iteration
%     end
%     X = Xn;
% 
% end
% 
% X_opt = X;
% F0 = sum(beta.*gamma.*X)';
% W0 = zeros(K,1);
% for k = 1:K
%     cosnt=0;
%     for ii=1:K
%         cosnt = cosnt + C(k,ii)*sum(X(:,ii).*beta(:,k).*gamma(:,ii));
%     end
%     W0(k)=1 + T_u*Pd*cosnt;
% end
% zsqr0 = sum(gamma.*sqrt(X))';
% Z0 = zsqr0.^2 + F0;
% A0 = F0 ./ (W0 .* Z0);
% t_opt = max(A0); %Compute t with optimized X
% 
% Maxnor_Var_pgd(iii) = t_opt; 
% 
% Pr0 = ((noise_p)* T_u*Pd* (X_opt.*Gammaa));
% for q = 1:K
%     b_pgd(q) = b_num(:,q)'*Pr0(:,q);
%     A_pgd(q) = Pr0(:,q)'*A_den(:,:,q)*Pr0(:,q);
%     Sig_pgd(q) = b_pgd(q)/A_pgd(q);
% end
% L_mse0_pgd = Sig_pgd';
% L_mse_ave_pgd(iii) = ( L_mse0_pgd(1) + L_mse0_pgd(2) ) /2; 
% 
% Pr_pgd(iii) = sum(sum(Pr0)); 
% 
% 
% for k = 1:K
%         for kk = 1:K
%             Skkk_pgd(k,kk) = X_opt(:,kk)'*(BETAA(:,k).*Gammaa(:,kk)); %Skk'
%         end
%     end
%  for k = 1:K
%         Var_minmax_pgd(k) = sum(X_opt(:,k).*BETAA(:,k).*Gammaa(:,k)) ./ (1+T_u*Pd* (Skkk_pgd(k,:)*PhiPhi2(k,:)') ); %min-max variance for each k
%  end
% 
% 
% for k = 1:K %Order of the 'correct' UE
%         for ii = 1:K %Order of all UE
%             a_pgd(k,ii) = X_opt(:,ii)'*(BETAA(:,k).*Gammaa(:,ii)) + PhiPhi2(k,ii).*( ( X_opt(:,ii).^(1/2) )'*(BETAA(:,k)./BETAA(:,ii).*Gammaa(:,ii)) ).^2; %E{|channel|^2}
%         end
%     end
%     for k = 1:K
%         nor_pgd(k) = Var_minmax_pgd(k)/a_pgd(k,k);
%     end
% 
%  [max_pgd index_pgd] = max(nor_pgd); 
% 
% %Zero-forcing SNR calculation - PGD
% P_u_pgd = noise_p*Pd* M; %Power to the UE 
% BETAA_u = sum(BETAA); 
% for k = 1:K %Variance of estimated channel
%     Kap_pgd(k) = ( T_u*Pd* ( Skkk_pgd(k,k) )^2 ) / (1+T_u*Pd* (Skkk_pgd(k,:)*PhiPhi2(k,:)') );
% end
% for k = 1:K %E{|w_ZF|^2}*Power
%     xsi_pgd(k) = ( BETAA_u(k)/(  Kap_pgd(k)*  (M-2) )) * P_u_pgd; 
% end
% for k = 1:K
%     for ii = 1:K
%         SNR_pgd(k,ii) = xsi_pgd(k)/( xsi_pgd(ii)/BETAA_u(ii)*Var_minmax_pgd(k) +    noise_p   );
%     end
% end
% %Average for two UE
% SNR_ave_pgd(iii) = (SNR_pgd(1,2) + SNR_pgd(2,1))/2; 
% 
% 
% if index_pgd == 1
%     SNR_ue_pgd(iii) = SNR_pgd(1,2);
% else
%     SNR_ue_pgd(iii) = SNR_pgd(2,1);
% end
% 
% 
% Var_ini(iii) = Var_mk3; %The worst channel estimation without DL

%% Pilot Asignment: (random choice) for extended UL traning
[U1,S,V]=svd(randn(TAUu,TAUu));%U includes tau orthogonal sequences
Phii=zeros(TAUu,K);

Point=randperm(TAUu,K); %If the number of UE is less than traning sequences, then no repeat
for k=1:K  
    Phii(:,k)=U1(:,Point(k));
end

Phii_cf333 = Phii; % pilot set of cell-free systems

% Create Gamma matrix for extended UL traning (variances of the channel estimates)
Gammaa333 = zeros(M,K);
ma333u=zeros(M,K);
for m=1:M
    for k=1:K
        mau333(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf333(:,k)'*Phii_cf333))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa333(m,k)=TAUu*Pu*BETAA(m,k)^2/(TAUu*Pu*mau333(m,k) + 1);
    end
end

Var_mk00333 = (BETAA - Gammaa333)./BETAA;
Var_mk0333 = sum(Var_mk00333) ./ M;
[Var_mk333 index_mk333] = max(Var_mk0333);

Var_mk333_ini(iii) = Var_mk333;


end

Maxnor_Var_iii(z) = sum(Maxnor_Var)/ii_oo; %average of iii times at this POSITION
%Maxnor_Var_pgd_iii(z) = sum(Maxnor_Var_pgd)/ii_oo; 
Maxnor_full_iii(z) = sum(Maxnor_full)/ii_oo; 

L_mse_ave_iii(z) = sum(L_mse_ave)/ii_oo; 
%L_mse_ave_pgd_iii(z) = sum(L_mse_ave_pgd)/ii_oo; 
L_mse_full_iii(z) = sum(L_mse_full)/ii_oo; 
%L_mse000_iii(z) = sum(L_mse000)/ii_oo; 

%Var_ini_iii(z) = Var_ini/ii_oo; 
Var_mk333_iii(z) = Var_mk333_ini/ii_oo; 


Pr00_iii(z) = sum(Pr00)/ii_oo; 
Pr_cvx_iii(z) = sum(Pr_cvx)/ii_oo; 
%Pr_pgd_iii(z) = sum(Pr_pgd)/ii_oo; 
%Pr_rrr_iii(z) = sum(Pr_rrr)/ii_oo; 

SNR_full_iii(z) = sum(SNR_ave_full)/ii_oo; 
SNR_sca_iii(z) = sum(SNR_ave_sca)/ii_oo; 
%SNR_pgd_iii(z) = sum(SNR_ave_pgd)/ii_oo; 
SNR_wDL_iii(z) = sum(SNR_ave_wDL)/ii_oo; 

SNR_ue_full_iii(z) = sum(SNR_ue_full)/ii_oo;
SNR_ue_sca_iii(z) = sum(SNR_ue_sca)/ii_oo;
%SNR_ue_pgd_iii(z) = sum(SNR_ue_pgd)/ii_oo;
SNR_ue_wDL_iii(z) = sum(SNR_ue_wDL)/ii_oo;

end
Maxnor_full_z(nn,:) = Maxnor_full_iii; 
L_mse_full_z(nn,:) = L_mse_full_iii; 

 Maxnor_Var_z(nn,:) = Maxnor_Var_iii; 
 L_mse_ave_z(nn,:) = L_mse_ave_iii; 
% 
% Maxnor_Var_pgd_z(nn,:) = Maxnor_Var_pgd_iii; 
% L_mse_ave_pgd_z(nn,:) = L_mse_ave_pgd_iii; 

%L_mse000_z(nn,:) = L_mse000_iii; 

%Var_ini_z(nn,:) = Var_ini_iii; 
Var_mk333_z(nn,:) = Var_mk333_iii; 

Pr00_z(nn,:) = Pr00_iii; 
Pr_cvx_z(nn,:) = Pr_cvx_iii; 
%Pr_pgd_z(nn,:) = Pr_pgd_iii; 
%Pr_rrr_z(nn,:) = Pr_rrr_iii; 


SNR_full_z(nn,:) = SNR_full_iii; 
SNR_sca_z(nn,:) = SNR_sca_iii; 
%SNR_pgd_z(nn,:) = SNR_pgd_iii; 
SNR_wDL_z(nn,:) = SNR_wDL_iii; 

SNR_ue_full_z(nn,:) = SNR_ue_full_iii;
SNR_ue_sca_z(nn,:) = SNR_ue_sca_iii;
%SNR_ue_pgd_z(nn,:) = SNR_ue_pgd_iii;
SNR_ue_wDL_z(nn,:) = SNR_ue_wDL_iii;

end

Maxnor_full_nn = sum(Maxnor_full_z)/N; %Fullpower
L_mse_full_nn = sum(L_mse_full_z)/N; 

Maxnor_Var_nn = sum(Maxnor_Var_z)/N; %CVX
L_mse_ave_nn = sum(L_mse_ave_z)/N; 
% 
% Maxnor_Var_pgd_nn = sum(Maxnor_Var_pgd_z)/N; %PGD
% L_mse_ave_pgd_nn = sum(L_mse_ave_pgd_z)/N; 

%L_mse000_nn = sum(L_mse000_z)/N; %Same radr SNR

%Var_ini_nn = sum(Var_ini_z)/N; %Withot
Var_mk333_nn = sum(Var_mk333_z)/N; %Without extend


Pr00_nn = sum(Pr00_z)/N; 
Pr_cvx_nn = sum(Pr_cvx_z)/N; 
% Pr_pgd_nn = sum(Pr_pgd_z)/N; 
%Pr_rrr_nn = sum(Pr_rrr_z)/N; 

SNR_full_nn = sum(SNR_full_z)/N;
SNR_sca_nn = sum(SNR_sca_z)/N;
% SNR_pgd_nn = sum(SNR_pgd_z)/N;
SNR_wDL_nn = sum(SNR_wDL_z)/N;

SNR_ue_full_nn = sum(SNR_ue_full_z)/N;
SNR_ue_sca_nn = sum(SNR_ue_sca_z)/N;
% SNR_ue_pgd_nn = sum(SNR_ue_pgd_z)/N;
SNR_ue_wDL_nn = sum(SNR_ue_wDL_z)/N;

for i = 1:z


SNR_full(i) = ( (1 - Maxnor_full_nn(i) ) *Pd/(Pd+1) ) ;

%SCA
SNR_sca(i) = ( (1 - Maxnor_Var_nn(i) ) *Pd/(Pd+1) ) ; 

%PGD
%SNR_pgd(i) = ( (1 - Maxnor_Var_pgd_nn(i) ) *Pd/(Pd+1) ) ;

%wDL
SNR_wDL(i) = ( (1 - Var_mk333_nn(i) ) *Pd/(Pd+1) ) ;

end

save ppSNR.mat

