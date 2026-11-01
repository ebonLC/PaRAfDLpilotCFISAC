clc
clear all

load('UE positions and results.mat')

%tau = 15; 
%T_u = 9; 

    Maxnor_full = [ ]; 

%% Selected Orthogonalties
for repst = 1:10

inde = [15 12 10 8 6 5 4 3 2 1];     

%% Initialization

%Colerration sequences generation
rho = 10^(-inde(repst)); % desired pairwise inner‐product
R = rho * ones(K,K) + (1-rho)*eye(K);
L = chol(R,'lower'); % L is K×K lower‐triangular, R = L * L'
X0 = L'; % size K×K, each column unit‐norm, X0' * X0 = R
Xpad = [ X0 ; zeros(tau-K, K) ];   % size tau x K, but rows K+1:10 are all zeros( R^N by padding zeros below
A = randn(tau); % random matrix, i.i.d. Gaussian
[U,~] = qr(A,0); % “economy” QR: U is orthonormal
X = U * Xpad; %Left‐multiply Xpad by U to “rotate” into a fully‐populated 10×5 matrix

Phii_cf = X; % pilot set of cell-free systems
Phii_ULr = Phii_cf; %Record, used for conventional

% Create Gamma matrix (variances of the channel estimates)
Gammaa = zeros(M,K);
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)= (norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf)   ))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau*Pu*BETAA(m,k)^2/(tau*Pu*mau(m,k) + 1);
    end
end

PhiPhi=zeros(K,K);
for k=1:K
    for ii=1:K
        PhiPhi(k,ii)=norm(Phii_cf(:,k)'*Phii_cf(:,ii));
        PhiPhi2(k,ii)=PhiPhi(k,ii)^2;
    end
end

%When without DL, the UL estimation error is the error for the whole system
Var_mk003 = (BETAA - Gammaa) ./ BETAA; %normalized
Var_mk03 = sum(Var_mk003) ./ M; %system average to each UE
Var_mk3 = max(Var_mk03); %normalized worst average value

Var_mk3_w(repst) = Var_mk3;

Var_real_w(repst) = max( sum(BETAA - Gammaa) ./ M) ;

%% Pilot Asignment: (random choice) for DL

X = []; R = []; L = []; X0 = []; Xpad = []; A = []; U = []; 
%Colerration sequences generation
rho = 10^(-inde(repst)); % desired pairwise inner‐product
R = rho * ones(K,K) + (1-rho)*eye(K);
L = chol(R,'lower'); % L is K×K lower‐triangular, R = L * L'
X0 = L'; % size K×K, each column unit‐norm, X0' * X0 = R
Xpad = [ X0 ; zeros(T_u-K, K) ];   % size tau x K, but rows K+1:10 are all zeros( R^N by padding zeros below
A = randn(T_u); % random matrix, i.i.d. Gaussian
[U,~] = qr(A,0); % “economy” QR: U is orthonormal
X = U * Xpad; %Left‐multiply Xpad by U to “rotate” into a fully‐populated 10×5 matrix

Phii_cf12 = X; % pilot set of cell-free systems

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
L_mse_full(repst) = sum(L_mse_f) /K; 

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
   Un_Var(k) =  T_u*Pd* ([Skk_full(k,[1:k-1]) Skk_full(k,[k+1:K])]*  [PhiPhi2(k,[1:k-1]) PhiPhi2(k,[k+1:K])]' ); 
end
    for k = 1:K
        Var_minmax_full_good(k) = sum(Etaa(:,1).*BETAA(:,k).*Gammaa(:,k)) ./ (1+T_u*Pd* (Skk_full(k,:)*PhiPhi2(k,:)') ); %min-max variance for each k within good orhothogal
    end
    for k = 1:K
        Var_minmax_full(k) = ( sum(Etaa(:,1).*BETAA(:,k).*Gammaa(:,k)) + Un_Var(k) ) ./ (1+T_u*Pd* (Skk_full(k,:)*PhiPhi2(k,:)') ); %min-max variance for each k
    end
    for k = 1:K
        nor_full(k) = Var_minmax_full(k)/a_full(k,k);
    end

Maxnor_full(repst) = max(nor_full);

Max_full(repst) = max(Var_minmax_full);

%end

% x = 1:1:10; 
% 
% figure(1)
% semilogy(x,Var_real_w,'-o',x,Max_full,'-x')
% xticklabels({'10^{-15}','10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
% xlabel('Sequence average correlation coefficient')
% ylabel('Real Error')
% figure(2)
% semilogy(x,Var_mk3_w,'-o',x,Maxnor_full,'-x')
% xticklabels({'10^{-15}','10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
% xlabel('Sequence average correlation coefficient')
% ylabel('Normalized Error')

%% CVX within SCA+bisection

Xopt = []; %initialize Xopt, for judging if the CVX is finally feasible
Xopt0 = [];

for m=1:M
    Xopt0(m)=1/(sum(Gammaa(m,:)));
end
Xopt = [Xopt0' Xopt0' Xopt0' Xopt0' Xopt0'];
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
    Maxnor_Var(repst) = max(nor_Var); 
else %If problem has solution, just tnext is the answer, thenjust calculate Pr
    Pr = ((noise_p)* T_u*Pd* (Xopt.*Gammaa));
    for k = 1:K %Order of the 'correct' UE
        for ii = 1:K %Order of all UE
            a_kr(k,ii) = Xopt(:,ii)'*(BETAA(:,k).*Gammaa(:,ii)) + PhiPhi2(k,ii).*( ( Xopt(:,ii).^(1/2) )'*(BETAA(:,k)./BETAA(:,ii).*Gammaa(:,ii)) ).^2; %E{|channel|^2}
        end
    end
    Maxnor_Var(repst) = tnext; 
end

for q = 1:K
    b(q) = b_num(:,q)'*Pr(:,q);
    A(q) = Pr(:,q)'*A_den(:,:,q)*Pr(:,q);
    Sig(q) = b(q)/A(q);
end
L_mse0 = Sig';
L_mse_ave(repst) = sum( L_mse0 ) /2; 

Pr_cvx(repst) = sum(sum(Pr)); 

%% PGD within feasible projection

%Avoid intereference
gamma = Gammaa; 
beta = BETAA; 
C = PhiPhi2; 
Desc   = 1e-15; %Gradient step size
maxite = 1000;
Fcc = noise_p * T_u * Pd;

X = ones(M,K) * (1/mean(gamma(:))); %Reset X
X = max(X,0);
for m = 1:M
    s = gamma(m,:)*X(m,:)';
    if s > 1
        X(m,:) = X(m,:) * (1/s); %(28c)
    end
end

for it = 1:maxite
    %F, W, Z as CVX
    F = sum(beta.*gamma.*X)';  
    W = zeros(K,1);
    for k = 1:K
        cosnt = 0;
        for ii = 1:K
            cosnt = cosnt + C(k,ii)*sum(X(:,ii).*beta(:,k).*gamma(:,ii));
        end
        W(k) = 1 + T_u*Pd * cosnt;
    end
    zsqr = sum(gamma.*sqrt(X),1)';
    Z = zsqr.^2 + F;
    A = F ./ (W .* Z);               

    %Worst A, the same as derivation, but notice gradiant A
    [~, kmax] = max(A);
    dF    = beta(:,kmax).*gamma(:,kmax);
    dzsqr    = gamma(:,kmax)./(2*sqrt(X(:,kmax)   +  1e-12 ));
    dZ    = 2*zsqr(kmax)*dzsqr + dF;
    gradA = zeros(M,K);
    gradA(:,kmax) = ( dF*Z(kmax) - F(kmax)*dZ ) / ( W(kmax)*Z(kmax)^2 ); %Quotient rule, dy/dx = ..../v^2

    %Descent step
    Xn = X - Desc*gradA; 
    Xn = real(Xn);
    Xn = max(Xn, 0);

    %Constraint/ feasible projection
for Rroww = 1: 15 %Xn on each projection, try to get the Xn close to three at the same time    
    Xn = max(Xn, 0); %(28b)

    for p = 1:K
        y = Xn(:,p).*gamma(:,p);
        for i = 1:M
            a     = thre * Fcc * A_den(i,:,p)';
            viol  = b_num(i,p) - a'*y;
            if viol > 0
                y = y + (viol/(a'*a)) * a; %(28d)
            end
        end
        Xn(:,p) = y ./ gamma(:,p);
    end

    for m = 1:M
        u = Xn(m,:)';  g = gamma(m,:)';
        s = g'*u;
        if s > 1
            u = u - ((s-1)/(g'*g))*g; %(28c), halfspace
            
            u_s = min(u);
            if u_s < 0
                for m = 1:M
                    s_s = gamma(m,:)*Xn(m,:)';
                    if s_s > 1
                        Xn(m,:) = Xn(m,:) * (1/s_s); %Avoid minus X and over power
                        u = Xn(m,:)'; 
                    end
                end
            end
        end
        Xn(m,:) = u';
    end
end
%If the iteration still cannot avoid negative or complex X...
Xn = real(Xn);
Xn = max(Xn,0); 

    if norm(Xn - X,'fro')/norm(X,'fro') < epsi/100
        X = Xn;
        break; %No need to reach the full iteration
    end
    X = Xn;

end

X_opt = X;
F0 = sum(beta.*gamma.*X)';
W0 = zeros(K,1);
for k = 1:K
    cosnt=0;
    for ii=1:K
        cosnt = cosnt + C(k,ii)*sum(X(:,ii).*beta(:,k).*gamma(:,ii));
    end
    W0(k)=1 + T_u*Pd*cosnt;
end
zsqr0 = sum(gamma.*sqrt(X))';
Z0 = zsqr0.^2 + F0;
A0 = F0 ./ (W0 .* Z0);
t_opt = max(A0); %Compute t with optimized X

Maxnor_Var_pgd(repst) = t_opt; 

Pr0 = ((noise_p)* T_u*Pd* (X_opt.*Gammaa));
for q = 1:K
    b_pgd(q) = b_num(:,q)'*Pr0(:,q);
    A_pgd(q) = Pr0(:,q)'*A_den(:,:,q)*Pr0(:,q);
    Sig_pgd(q) = b_pgd(q)/A_pgd(q);
end
L_mse0_pgd = Sig_pgd';
L_mse_ave_pgd(repst) = ( L_mse0_pgd(1) + L_mse0_pgd(2) ) /2; 

Pr_pgd(repst) = sum(sum(Pr0)); 

Var_ini(repst) = Var_mk3; %The worst channel estimation without DL

%% Pilot Asignment: (random choice) for extended UL traning
TAUu = tau + T_u; 
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
Var_mk333 = max(Var_mk0333);

Var_mk333_ini(repst) = Var_mk333;



end


