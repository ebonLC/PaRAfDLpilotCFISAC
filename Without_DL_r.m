function [Xopt,L_mse00,KK,Gammaaa,Pr_rr] = Without_DL_r(SNR_aver,M,K,Phii_ULr,tau,Pu,Pd,BETAA,Gammaa,b_num,A_den,noise_p,alpha,h_los,T_u)

L1 = 1;
Phii_cf = Phii_ULr;

%% 1) Each AP has equal power allocations for K terminals
%%Compute etaa(m): (each AP transmits equal power to K terminals)
etaa=zeros(M,1);
for m=1:M
    etaa(m)=1/(sum(Gammaa(m,:)));
end

%%Compute Rate
SINR=zeros(1,K);
R_cf=zeros(1,K);

%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k)=sum((etaa.^(1/2)).*((Gammaa(:,ii)./BETAA(:,ii)).*BETAA(:,k)))*Phii_cf(:,ii)'*Phii_cf(:,k);
    end
end
PC1=(abs(PC)).^2;

for k=1:K
    num=0;
    for m=1:M
        num=num + (etaa(m)^(1/2))*Gammaa(m,k);
    end
    SINR(k) = Pd*num^2/(1 + Pd*sum(BETAA(:,k)) + Pd*sum(PC1(:,k)) - Pd*PC1(k,k));
    %Rate:
    R_cf(k) = log2(1+ SINR(k));
end

R_cf_min = min(R_cf);

%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin= max(2^(R_cf_min-2)-1 , 0);
tmax=2^(2*R_cf_min+2.5)-1;
epsi=max(tmin/5,0.001);
KK = K+K;
PhiPhi=zeros(K,K);
X_r = ones(M,K);
for ii=1:K
    for k=1:K
        PhiPhi(ii,k)=norm(Phii_cf(:,ii)'*Phii_cf(:,k));
    end
end
BETAAn=BETAA*Pd;
Gammaaa = [Gammaa Gammaa];
Gammaan=Gammaaa*Pd;

for m=1:M
    Xopt0(m)=1/(sum(Gammaaa(m,:))); %assume each AP has equal power to UE
end
Xopt = [sqrt(Xopt0') sqrt(Xopt0') sqrt(Xopt0') sqrt(Xopt0')];
yopt=sqrt(Pd)*ones(M,1);
Zopt=ones(K,K);
y0opt=sqrt(Pd)*ones(M,1);
Z0opt=ones(K,K);
%cvx_solver sedumi
count = 0;
cON_r = 1;
Cou_r = 3;
cvx_quiet true
    while cON_r <= Cou_r
        cON_r
            while( tmax - tmin > epsi)
            count = count + 1;
            tnext = (tmax+tmin)/2;
           cvx_begin sdp
              variables X(M,KK) y(M,1) Z(K,K) y0(M,1) Z0(K,K) SNRr(K,1) %%Notice here X is sqrt(eta)
              expression SNR_Taylor(M,M,K)
              minimize(0)
              subject to
                for k=1:K
                     norm([sqrt(tnext)*[Z((1:(k-1)),k);Z(((k+1):K),k)].*[PhiPhi((1:(k-1)),k);PhiPhi(((k+1):K),k)];... %user-interference
                         sqrt(tnext/L1)*(BETAAn(:,k).^(1/2)).*y ; ... %estimation uncertainty
                         sqrt(tnext)*[Z0(:,k)].*[PhiPhi(:,k)]; ... %sensing-interference PI
                         sqrt(tnext/L1)*(BETAAn(:,k).^(1/2)).*y0 ; ... %sensing-interference PII
                         sqrt(tnext*Pd)/L1]) <= (Gammaan(:,k))'*X(:,k) ;
                end
                for m=1:M
                    norm(((Gammaan(m,1:K)).^(1/2)).*X(m,1:K)) <= y(m); %BU, channel estimation uncertainty
                    %y(m)<=(sqrt(Pd)/L1);
                end                
                for k=1:K
                    for ii=1:K
                        sum( ((Gammaan(:,ii)./BETAA(:,ii)).*BETAA(:,k)).*X(:,ii)  ) <= Z(ii,k); %UI, other-UE interference, notice that "k'≠k"eta*gamma*beta parts are moved to with BU
                    end
                end
                %Caused by sensing interference
                for m=1:M
                    norm(((Gammaan(m,1:K)).^(1/2)).*X(m,(K+1):KK)) <= y0(m);
                    %y0(m)<=(sqrt(Pd)/L1);
                end
                for k=1:K
                    for ii=1:K
                        sum( ((Gammaan(:,ii)./BETAA(:,ii)).*BETAA(:,k)).*X(:,(K+ii))  ) <= Z0(ii,k);
                    end
                end
                Pr = ((noise_p))* T_u*Pd* ( (X(:,(K+1:KK)).^2) .*Gammaaa(:,(K+1:KK)) ) ; %Original
                Pr_sqt = (((noise_p))* T_u*Pd).^(1/2) *  ( X(:,(K+1:KK)) .*(Gammaaa(:,(K+1:KK)).^(1/2)) ); %sqrt() expression, which is the x
                Pr_sqt_opt = (((noise_p))* T_u*Pd).^(1/2) *  ( Xopt(:,(K+1:KK)) .*(Gammaaa(:,(K+1:KK)).^(1/2)) ); %x0 expression, which is the input x0
                for i = 1:M %Transmit AP
                    for j = 1:M %Receive AP
                        for k = 1:K %User
                            SNR_Taylor(i,j,k) = (2*alpha(i,j,k).* (h_los(i,j)) .*Pr_sqt_opt(i,k).*Pr_sqt(i,k) - alpha(i,j,k).* (h_los(i,j)) .*(Pr_sqt_opt(i,k).^2)) / (noise_p);
                            %SNR_Taylor(i,j,k) >= 1e-9; %Constraints,make the radar SINR not less than the compared model
                        end
                    end
                end
                reshape( sum(sum(SNR_Taylor)), [2,1] ) == max(SNR_aver); %Constraints,make the radar SINR not less than the compared model
                 for m=1:M
                    for k=1:KK
                        X(m,k)>=0;
                    end
                end
                for m=1:M
                    norm(((Gammaan(m,:)).^(1/2)).*X(m,:)) <= (sqrt(Pd)/L1); %Power constraint for each AP
                end
            cvx_end
            
            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            t_r = tmin;
            Xopt=X; %Also let Pr_sqt_opt = Pr_sqt
            yopt=y;
            Zopt=Z;
            y0opt = y0;
            Z0opt = Z0;
            SNR_Taylor_r = SNR_Taylor; %Record SNR
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end
        end
            cON_r = cON_r+1;
            tmin= max(2^(R_cf_min-2)-1 , 0);
            tmax=2^(2*R_cf_min+2.5)-1;
    end

%R_cf_opt = log2(1+t_r);

Pr_r = (noise_p)* T_u*Pd* ( (Xopt(:,(K+1:KK)).^2) .*Gammaaa(:,(K+1:KK)) );
for q = 1:K
    b(q) = b_num(:,q)'*Pr_r(:,q);
    A(q) = Pr_r(:,q)'*A_den(:,:,q)*Pr_r(:,q);
    Sig(q) = b(q)/A(q);
end
L_mse00 = Sig';
Pr_rr = (noise_p)* Pd* ( (Xopt(:,(K+1:KK)).^2) .*Gammaaa(:,(K+1:KK)) );
count