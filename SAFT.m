%% SAFT-gamma Mie Code
%Directory:
%>>SAFT.m (l.73)
    %Global parameters calculated:
        %>>SegmentDensity.m (l.178): Finds the segment density (rho_s)
        %>>MieC.m (l.190): Coefficient of the Mie potential (C)
        %>>HSdiameter.m (l.200): Calculates the reference hard sphere diameter for
        %each combination of group (HSd)
        %>>SegmentFraction.m (l.211): Finds the segment fractions (xSegment)
        % Calculated in zetaFunction.m
        %>>x0.m (l.232): Finds the ratios of sigma to HSd for each group pair (x_0)
        %Calculated in A1.m
        %>>zetaStar.m (l.242): Finds the packing fraction (zeta2) 
        %Calculated in CorrFrac.m
        %>>GroupFracComp.m (l.253): Finds the group fraction on a component (z)
        %Calculated in Chain.m
        %>>EffectiveDiameter.m (l.267): Finds the effective diameter cubed (EffD3)
        %Calculated in Chain.m
        %>>EffectiveSigma.m (l.284): Calculates the effective sigma (EffSigma)
        %Calculated in Chain.m
        %>>zetaChain.m (l.301): Finds the packing fraction of the chain term 
        %(zetaC) Calculated in Chain.m
        %>>EffectiveC.m (l.323): Finds the effective MieC (EffC) Calculated in PDF1.m
    %>>Ideal.m (l.331): Finds the Ideal Term (no correlation) (A_I)
    %>>Ideal2.m (l.339): Finds the Ideal term using the correlation (A_I)
    %>>Mono.m (l.355): Mono Term (A_M)
        %>>HS.m (l.364): Hard-sphere contribution (A_HS)
            %>>zetaFunction.m (l.369): Finds the 4 moment densities (zdtaF)
        %>>A1.m (l.403): First perturbation (A_1)
            %>>Aa1.m (l.421): Finds the total specific A1 (a_1)
                %>>Aa1k.m (l.436): Finds the specific A1 for each group pair (a_1k)
                    %>>Aa1s.m (l.462): Finds the values a_1s for each lambda
                    %(a_1sA,a1sR)
                    %>>Bk.m (l.478): Finds the value of B_k for each lambda 
                    %(B_kA,B_kR)
        %>>A2.m (l.490): Second perturbation (A_2)
            %>>Aa2.m (l.504): Finds the total specific A2 (a_2)
                %>>Aa2k.m (l.517): Finds the specific A2 for each group pair (a_2k)
                    %>>CorrFrac.m (l.534): Correction Fraction (X)
                        %>>Falpha.m (l.547): Finds Fm for each alpha (F1,F2,F3)
        %>>A3.m (l.559): Third perturbation (A_3)
            %>>Aa3.m (l.573): Finds the total specific A3 (a_3)
                %>>Aa3k.m (l.583): Finds the specific A3 for each group pair (a_3k)
    %>>Assoc.m (l.597): Association Term (A_A)
        %>>AssociationFraction.m (l.620): Finds the fraction of association sites 
        %of a group on a component that are bonded (X_Assoc)
            %>>AssocStrength.m (l.646): Finds the 'strength' of a bond between
            %sites (Delta)
    %>>Chain.m (l.687): Chain Term (A_C)
        %>>EffectiveLambda.m (l.715): Calculates the effective lambdas (EffLambdaA,
        %(EffLambdaA,EffLambdaR)
        %>>EffectiveEpsilon.m (l.732): Calculates the effective epsilon
        %(EffEpsilon)
        %>>PDFMie.m (l.749): Finds the PDF for a mie fluid (gMie)
            %>>PDFHS.m (l.760): Finds the PDF for a HS fluid (gHS)
            %>>PDF1.m (l.774): Finds the first order PDF (g1)
                %>>EffectiveA1s.m (l.795): Finds the effective a1s (Effa1sA,Effa1sR)
                %>>EffectiveBk.m (l.808): Finds the effective Bk (EffBkA,EffBkR)
                %>>PartialA1.m (l.818): Finds the derivative of a1 with respect to
                %rho_s (Da1Drho_s)
                    %>>PartialA1s.m (l.829): Finds the derivative of a1s with
                    %respect to rho_s (Da1sDrho_sA,Da1sRDrho_sR)
                    %>>PartialBk.m (l.842): Finds the derivative of Bk with respect
                    %to rho_s (DBkDrho_sA,DBkDrho_sR)
            %>>PDF2.m (l.852): Finds the second order PDF (g2)
                %>>EffFac.m (l.863): Finds the Correction Factor (yC)
                %>>PDFMCA.m (l.872): Finds the second order macroscopic
                %compressibility approximation (gMCA)
                    %>>PartialA2.m (l.887): Finds the derivative of a2/(X-1) with
                    %respect to rho_s (Da2Drho_s)
%% Main Function
function [A_T]        = SAFT(NParticles,Vol,Temp,xComp,Components)  
global NComponent NGroup k_B h N_A rho_s C HSd R NSites
load Groups.mat
load CrossInteraction.mat
load NAssoc.mat
load Association_Epsilon.mat
load Bonding_Vol.mat

%Setting the important constants
h=6.62607004e-34;k_B=1.38064852e-23;N_A=6.02214086e23;R=8.314;

%Parameter Extraction
NComponent    = size(Components,2);
x             = 1;
for i=1:NComponent
    for n=1:13
        if Components(n,i)~=0
            G(x) = n;
            x    = x+1;
            n    = n+1;
        else
            n    = n+1;
        end
    end
end
H          = unique(G);
n_Assoc    = NAssoc(:,[H]);
Parameters = Groups(:,[H]);
NGroup     = size(Parameters,2);
Component  = Components([H],:);
EPSILON    = CrossInteraction([H],[H]);
y          = 1;

        
for i=1:NComponent
    for k=1:NGroup
        GMultiplicity(k,i) = (Component(k,i));
        M(k,i)             = (Parameters(1,k))*GMultiplicity(k,i);
        k                  = k+1;
    end
    Mr(i) = sum(M(:,i))/N_A;
    i     = i+1;
end
for k=1:NGroup
    lambdaR(k,k)     = (Parameters(3,k));
    lambdaA(k,k)     = (Parameters(2,k));
    sigma(k,k)       = (Parameters(4,k));
    NSegment(k)      = (Parameters(5,k));
    ShapeFactor(k)   = (Parameters(6,k));
    a(k)             = Parameters(11,k);
    b(k)             = Parameters(12,k);
    c(k)             = Parameters(13,k);
    d(k)             = Parameters(14,k);
    H0(k)            = Parameters(15,k);
    S0(k)            = Parameters(16,k);
    k                = k+1;
end
for k=1:NGroup
    for l=1:NGroup
        if k==l
            epsilon(k,l) = (EPSILON(k,l))*k_B;
        else
            sigma(k,l)   = (sigma(k,k)+sigma(l,l))/2;
            lambdaA(k,l) = 3+sqrt((lambdaA(k,k)-3)*(lambdaA(l,l)-3));
            lambdaR(k,l) = 3+sqrt((lambdaR(k,k)-3)*(lambdaR(l,l)-3));
            if EPSILON(k,l)==0
                epsilon(k,l) = sqrt(sigma(k,k)^3*sigma(l,l)^3)*sqrt(EPSILON(k,k)*EPSILON(l,l))/sigma(k,l)^3*k_B;
            else
                epsilon(k,l) = (EPSILON(k,l))*k_B;
            end
        end
    end
    k = k+1;
end 

NSites=3;
n_Assoc       = NAssoc(:,[H]');
BondVol       = Bonding_Vol([H],[H],:,:);
AssocEpsilon  = Association_Epsilon([H],[H],:,:)*k_B;
                

%Defining Global Variables
rho_s         = SegmentDensity(NParticles/Vol,GMultiplicity,NSegment,ShapeFactor,xComp);
C             = MieC(lambdaR,lambdaA);
HSd           = HSdiameter(lambdaR,lambdaA,sigma,epsilon,Temp);

%Calling each contribution
% A_I           = Ideal(NParticles,Vol,Temp,xComp,Mr);
A_I           = Ideal2(NParticles,Temp,Vol,a,b,c,d,H0,S0,xComp,GMultiplicity);
A_M           = Mono(xComp,lambdaR,lambdaA,epsilon,sigma,GMultiplicity,NSegment,ShapeFactor,NParticles,Temp,Vol);
A_C           = Chain(lambdaR,lambdaA,epsilon,sigma,GMultiplicity,NSegment,ShapeFactor,Temp,xComp,NParticles,Vol);
A_A           = Assoc(xComp,GMultiplicity,NParticles,n_Assoc,BondVol,AssocEpsilon,Temp,Vol,lambdaR,epsilon,sigma);
A_T           = A_I+A_M+A_A+A_C;

%Note: Regular A, B, C and D are used to find sums throughout
%Note: i is used to iterate through components, k and l are used to
%iterate through groups and a and b are used to iterate through association
%sites
%Note: This code uses SI units:
    %Vol: m^3, Pressure: Pa, Temp: K, Energy: J, Mass: kg
%Note: Diagnostic.m can be used to find the AAD to verify the code
%Note: Main_Pure.m and Main_Binary.m produce PV diagrams
%Note: Test.m can be used to Test individual functions in isolation
end
%% Global Variables
function [rho_s]      = SegmentDensity(rho,GMultiplicity,NSegment,ShapeFactor,xComp)
global NGroup NComponent
for n=1:NComponent
    for k=1:NGroup
        A(k) = GMultiplicity(k,n)*NSegment(k)*ShapeFactor(k);
        k    = k+1;
    end
    B(n) = xComp(n)*sum(A);
    n    = n+1;
end
rho_s = rho*sum(B);
end
function [C]          = MieC(lambdaR,lambdaA)
global NGroup
for k=1:NGroup
    for l=1:NGroup
        C(k,l) = (lambdaR(k,l)./(lambdaR(k,l)-lambdaA(k,l))).*(lambdaR(k,l)./lambdaA(k,l)).^(lambdaA(k,l)./(lambdaR(k,l)-lambdaA(k,l)));
        l      = l+1;
    end
    k = k+1;
end
end
function [d]          = HSdiameter(lambdaR,lambdaA,sigma,epsilon,Temp)
global C k_B NGroup
for k=1:NGroup
    for l=1:NGroup
        Mie    = @(r) 1-exp(-(C(k,l).*epsilon(k,l).*((sigma(k,l)./r).^lambdaR(k,l)-(sigma(k,l)./r).^lambdaA(k,l)))./(k_B*Temp));
        d(k,l) = integral(Mie,0,sigma(k,l));
        l      = l+1;
    end
    k = k+1;
end
end
function [xSegment]   = SegmentFraction(xComp,GMultiplicity,NSegment,ShapeFactor)
global NGroup NComponent
for i=1:NComponent
    for k=1:NGroup
        A(k) = GMultiplicity(k,i).*NSegment(k).*ShapeFactor(k);
        k    = k+1;
    end
    B(i) = xComp(i).*sum(A);
    i    = i+1;
end
C           = sum(B);
for k=1:NGroup
for i=1:NComponent
    D(i) = xComp(i).*GMultiplicity(k,i).*NSegment(k).*ShapeFactor(k);
    i    = i+1;
end

xSegment(k) = sum(D)/C;
k           = k+1;
end
end
function [x_0]        = x0(sig)
global HSd NGroup
for k=1:NGroup
    for l=1:NGroup
        x_0(k,l) = sig(k,l)/HSd(k,l);
        l        = l+1;
    end
    k = k+1;
end
end
function [zeta2]      = zetaStar(sigma)
global xSegment NGroup rho_s
for k=1:NGroup
    for l=1:NGroup
        A(k,l) = xSegment(k)*xSegment(l)*sigma(k,l)^3;
        l      = l+1;
    end
    k = k+1;
end
zeta2=rho_s*pi*sum(sum(A))/6;
end
function [zEff]       = GroupFracComp(GMultiplicity,NSegment,ShapeFactor)
global NComponent NGroup
for i=1:NComponent
    for k=1:NGroup
        A(k)  = GMultiplicity(k,i)*NSegment(k)*ShapeFactor(k);
        k     = k+1;
    end
    for l=1:NGroup
        zEff(l,i) = A(l)/sum(A);
        l         = l+1;
    end
    i = i+1;
end
end
function [EffD3]      = EffectiveDiameter()
global z HSd NGroup NComponent
for i=1:NComponent
    for j=1:NComponent
        for k=1:NGroup
            for l=1:NGroup
                A(k,l) = z(k,i)*z(l,j)*HSd(k,l)^3;
                l      = l+1;
            end
            k = k+1;
        end
        EffD3(i,j) = sum(sum(A));
        j          = j+1;
    end
    i = i+1;
end
end
function [EffSigma]   = EffectiveSigma(sigma)
global NComponent NGroup z
for i=1:NComponent
    for j=1:NComponent
        for k=1:NGroup
            for l=1:NGroup
                A(k,l) = z(k,i)*z(l,j)*sigma(k,l)^3;
                l      = l+1;
            end
            k = k+1;
        end
        EffSigma(i,j) = (sum(sum(A)))^(1/3);
        j              = j+1;
    end
    i = i+1;
end
end
function [zetaC]      = zetaChain(xComp,GMultiplicity,NSegment,ShapeFactor,rho_s1)
global NComponent NGroup EffD3 zeta1 rho_s

for i=1:NComponent
    for j=1:NComponent
        for k=1:NGroup
            A(k) = NSegment(k)*GMultiplicity(k,i)*ShapeFactor(k);
            k    = k+1;
        end
        for l=1:NGroup
            B(l) = NSegment(l)*GMultiplicity(l,j)*ShapeFactor(l);
            l    = l+1;
        end
        C(j) = xComp(i)*sum(A)*xComp(j)*sum(B)*EffD3(i,j);
        j    = j+1;
    end
    i = i+1;
end

zetaC = zeta1*rho_s1/rho_s;

end
function [EffC]       = EffectiveC(EffLambdaR,EffLambdaA)
global NComponent
for i=1:NComponent
    EffC(i,i) = (EffLambdaR(i,i)./(EffLambdaR(i,i)-EffLambdaA(i,i))).*(EffLambdaR(i,i)./EffLambdaA(i,i)).^(EffLambdaA(i,i)./(EffLambdaR(i,i)-EffLambdaA(i,i)));
    i         = i+1;
end
end
%% Ideal Term
function [A_I]        = Ideal(NParticles,Vol,Temp,xComp,Mr)
global h k_B NComponent
for i=1:NComponent
    A(i,1) = xComp(i)*log(xComp(i)*NParticles*(h^2/(2*Mr(i)*pi*k_B*Temp))^(3/2)/Vol);
    i      = i+1;
end
A_I = (sum(A)-1)*NParticles*k_B*Temp;
end
function [A_I]        = Ideal2(NParticles,Temp,Vol,a,b,c,d,H0,S0,xComp,GMultiplicity)
global N_A NGroup R NComponent
mole = NParticles/N_A;
for i=1:NComponent
    for k=1:NGroup
        Cp1  = @(T) a(k)+b(k)*T+c(k)*T.^2+d(k)*T.^3;
        Cp2  = @(T) a(k)./T+b(k)+c(k)*T+d(k)*T.^2;
        A(k) = (integral(Cp1,298,Temp)+H0(k)-R*Temp-Temp*integral(Cp2,298,Temp)-S0(k)*Temp-Temp*R*log(Vol/0.02477572))*GMultiplicity(k,i);
        k    = k+1;
    end
    B(i) = xComp(i)*sum(A);
    i    = i+1;
end
A_I = sum(B)*mole;
end
%% Monomer Term
function [A_M]        = Mono(xComp,lambdaR,lambdaA,epsilon,sigma,GMultiplicity,NSegment,ShapeFactor,NParticles,Temp,Vol)
global k_B
A_HS  = HS(xComp,sigma,GMultiplicity,NSegment,ShapeFactor,NParticles,Temp,Vol);
A_1   = A1(lambdaR,lambdaA,epsilon,sigma,NParticles,Vol,GMultiplicity,NSegment,ShapeFactor,xComp,Temp);
A_2   = A2(lambdaR,lambdaA,epsilon,sigma,NParticles,Temp,Vol,GMultiplicity,NSegment,ShapeFactor,xComp);
A_3   = A3(lambdaR,lambdaA,epsilon,NParticles,Temp,Vol,GMultiplicity,NSegment,ShapeFactor,xComp);
A_Dis = (A_1+A_2+A_3)*NParticles*Temp*k_B; %Quick fix with the perturbation terms (mimic paper exactly)
A_M   = A_HS*NParticles*Temp*k_B+A_Dis;
end
function [A_HS]       = HS(xComp,sig,GMultiplicity,NSegment,ShapeFactor,NParticles,Temp,Vol)
global k_B zetaF
zetaF                     = zetaFunction(xComp,GMultiplicity,NSegment,ShapeFactor);
A_HS                      = (6*Vol/(pi*NParticles))*(((zetaF(3).^3/zetaF(4).^2)-zetaF(1))*log(1-zetaF(4))+3*(zetaF(2).*zetaF(3))./(1-zetaF(4))+(zetaF(3).^3)/(zetaF(4).*(1-zetaF(4))^2));
end
function [zetaF]      = zetaFunction(xComp,GMultiplicity,NSegment,ShapeFactor)
global rho_s HSd NGroup xSegment 

xSegment = SegmentFraction(xComp,GMultiplicity,NSegment,ShapeFactor);

    for k=1:NGroup
        A(k) = xSegment(k).*HSd(k,k).^0;
        k    = k+1;
    end
    
    zetaF(1) = sum(A).*rho_s.*pi./6;
    
    for k=1:NGroup
        A(k) = xSegment(k).*HSd(k,k).^1;
        k    = k+1;
    end
    
    zetaF(2) =sum(A).*rho_s.*pi./6;
    
    for k=1:NGroup
        A(k) = xSegment(k).*HSd(k,k).^2;
        k    = k+1;
    end
    
    zetaF(3) = sum(A).*rho_s.*pi./6;
    
    for k=1:NGroup
        A(k) = xSegment(k).*HSd(k,k).^3;
        k    = k+1;
    end
    
    zetaF(4) = sum(A).*rho_s.*pi./6;
    
end
function [A_1]        = A1(lambdaR,lambdaA,epsilon,sigma,NParticles,Vol,GMultiplicity,NSegment,ShapeFactor,xComp,Temp)
global rho_s x_0 NComponent NGroup k_B
%Defining Global variable for A1 term
x_0 = x0(sigma);
for n=1:NComponent
    for k=1:NGroup
        A(k) = GMultiplicity(k,n)*NSegment(k)*ShapeFactor(k);
        k    = k+1;
    end
    B(n) = xComp(n)*sum(A);
    n    = n+1;
end
m_mean = sum(B);
%Calling the molecular helmholtz energy
a_1    = Aa1(lambdaR,lambdaA,epsilon);
%Finding the first perturbation term
A_1    = m_mean*a_1/(k_B*Temp);
end
function [a_1]        = Aa1(lambdaR,lambdaA,epsilon)
global NGroup xSegment

a_1k = Aa1k(lambdaR,lambdaA,epsilon);

for k=1:NGroup
    for l=1:NGroup
        A(l) = xSegment(k)*xSegment(l)*a_1k(k,l);
        l    = l+1;
    end
    B(k) = sum(A);
    k    = k+1;
end
a_1 = sum(B);
end
function [a_1k]       = Aa1k(lambdaR,lambdaA,epsilon)
global C NGroup x_0 xSegment HSd rho_s zeta1

for k=1:NGroup
    for l=1:NGroup
        A(l) = xSegment(k)*xSegment(l)*HSd(k,l)^3;
        l    = l+1;
    end
    B(k) = sum(A);
    k    = k+1;
end

zeta1 = pi*rho_s*sum(B)/6;
a_1sA = Aa1s(lambdaA,epsilon);
B_A   = Bk(lambdaA,epsilon);
a_1sR = Aa1s(lambdaR,epsilon);
B_R   = Bk(lambdaR,epsilon);

for k=1:NGroup
    for l=1:NGroup
        a_1k(k,l) = C(k,l)*(x_0(k,l)^(lambdaA(k,l))*(a_1sA(k,l)+B_A(k,l))-x_0(k,l)^(lambdaR(k,l))*(a_1sR(k,l)+B_R(k,l)));
        l         = l+1;
    end
    k = k+1;
end
end
function [a_1s]       = Aa1s(lambda,epsilon)
global rho_s NGroup HSd zeta1
for k=1:NGroup
    for l=1:NGroup
        A            = [ 0.81096   1.7888  -37.578   92.284;
                         1.02050 -19.3410  151.260 -463.500; 
                        -1.90570  22.8450 -228.140  973.920; 
                         1.08850  -6.1962  106.980 -677.640 ]*[1; 1/lambda(k,l); 1/lambda(k,l)^2; 1/lambda(k,l)^3];
        zeta_eff(k,l) =  A(1).*zeta1+A(2).*zeta1^2+A(3).*zeta1^3+A(4).*zeta1^4;
        a_1s(k,l)     =  -2*pi*rho_s*epsilon(k,l)*HSd(k,l)^3*(1-zeta_eff(k,l)/2)/((lambda(k,l)-3)*(1-zeta_eff(k,l))^3);
        a_1star(k,l)  =  a_1s(k,l)/epsilon(k,l);
        l             =  l+1;
    end
    k = k+1;
end
end
function [B_k]        = Bk(lambda,epsilon)
global x_0 zeta1 rho_s NGroup HSd
for k=1:NGroup
    for l=1:NGroup
        I        = (1-x_0(k,l).^(3-lambda(k,l)))/(lambda(k,l)-3);
        J        = (1-x_0(k,l).^(4-lambda(k,l))*(lambda(k,l)-3)+x_0(k,l)^(3-lambda(k,l))*(lambda(k,l)-4))/((lambda(k,l)-4)*(lambda(k,l)-3));
        B_k(k,l) = 2*pi*rho_s*HSd(k,l)^3*epsilon(k,l)*((1-zeta1/2)*I/(1-zeta1)^3-(9*zeta1*(1+zeta1)*J)/(2*(1-zeta1)^3));
        l        = l+1;
    end
    k = k+1;
end
end
function [A_2]        = A2(lambdaR,lambdaA,epsilon,sigma,NParticles,Temp,Vol,GMultiplicity,NSegment,ShapeFactor,xComp)
global rho_s k_B NComponent NGroup
a_2 = Aa2(lambdaR,lambdaA,epsilon,sigma);
for n=1:NComponent
    for k=1:NGroup
        A(k) = GMultiplicity(k,n)*NSegment(k)*ShapeFactor(k);
        k    = k+1;
    end
    B(n) = xComp(n)*sum(A);
    n    = n+1;
end
m_mean = sum(B);
A_2    = m_mean*a_2/(k_B*Temp)^2;
end
function [a_2]        = Aa2(lambdaR,lambdaA,epsilon,sigma)
global NGroup xSegment
a_2k = Aa2k(lambdaR,lambdaA,epsilon,sigma);
for k=1:NGroup
    for l=1:NGroup
        A(l) = xSegment(k)*xSegment(l)*a_2k(k,l);
        l    = l+1;
    end
    B(k) = sum(A);
    k    = k+1;
end
a_2  = sum(B);
end
function [a_2k]       = Aa2k(lambdaR,lambdaA,epsilon,sigma)
global zeta1 NGroup C x_0 K_HS
K_HS = (1-zeta1)^4/(1+4*zeta1+4*zeta1^2-4*zeta1^3+zeta1^4);
X    = CorrFrac(lambdaA,lambdaR,sigma);
        a_1s2A    = Aa1s(2*lambdaA,epsilon);
        a_1s2R    = Aa1s(2*lambdaR,epsilon);
        a_1sAR    = Aa1s(lambdaA+lambdaR,epsilon);
        
        B_k2A     = Bk(2*lambdaA,epsilon);
        B_k2R     = Bk(2*lambdaR,epsilon);
        B_kAR     = Bk(lambdaA+lambdaR,epsilon);
for k=1:NGroup
    for l=1:NGroup
        a_2k(k,l) = 0.5*K_HS*(1+X(k,l))*epsilon(k,l)*C(k,l)^2*(x_0(k,l)^(2*lambdaA(k,l))*(a_1s2A(k,l)+B_k2A(k,l))-2*x_0(k,l)^(lambdaR(k,l)+lambdaA(k,l))*(a_1sAR(k,l)+B_kAR(k,l))+x_0(k,l)^(2*lambdaR(k,l))*(a_1s2R(k,l)+B_k2R(k,l)));
    end
end
end
function [X]          = CorrFrac(lambdaA,lambdaR,sigma)
global rho_s NGroup C zeta2
zeta2 = zetaStar(sigma);
for k=1:NGroup
    for l=1:NGroup
        alpha  = C(k,l)*(1/(lambdaA(k,l)-3)-1/(lambdaR(k,l)-3));
        F1     = Falpha(alpha,1);
        F2     = Falpha(alpha,2);
        F3     = Falpha(alpha,3);
        X(k,l) = F1*zeta2+F2*zeta2^5+F3*zeta2^8;
    end
end      
end
function [Fm]         = Falpha(alpha,m)
Phi=[7.5365557,-359.4400,1550.9,-1.19932,-1911.28,9236.9,10;-37.60463,1825.6,-5070.1,9.063632,21390.175,-129430,10;71.745953,-3168,6534.6,-17.9482,-51320.7,357230,0.57;-46.83552,1884.2,-3288.7,11.34027,37064.54,-315530,-6.7;-2.467982,-0.82376,-2.7171,20.52142,1103.742,1390.2,-8;-0.50272,-3.1935,2.0883,-56.6377,-3264.61,-4518.2,NaN;8.0956883,3.709,0,40.53683,2556.181,4241.6,NaN];
for n=1:4
    A(n) = Phi(n,m)*alpha^(n-1);
    n    = n+1;
end
for j=5:7
    B(j) = Phi(j,m)*alpha^(j-4);
    j    = j+1;
end
Fm = sum(A)/(1+sum(B));
end
function [A_3]        = A3(lambdaR,lambdaA,epsilon,NParticles,Temp,Vol,GMultiplicity,NSegment,ShapeFactor,xComp)
global rho_s k_B NComponent NGroup
a_3 = Aa3(lambdaR,lambdaA,epsilon);
for n=1:NComponent
    for k=1:NGroup
        A(k) = GMultiplicity(k,n)*NSegment(k)*ShapeFactor(k);
        k    = k+1;
    end
    B(n) = xComp(n)*sum(A);
    n    = n+1;
end
m_mean = sum(B);
A_3    = m_mean*a_3/(k_B*Temp)^3;
end
function [a_3]        = Aa3(lambdaR,lambdaA,epsilon)
global NGroup xSegment
a_3k = Aa3k(lambdaR,lambdaA,epsilon);
for k=1:NGroup
    for l=1:NGroup
        A(k,l) = xSegment(k)*xSegment(l)*a_3k(k,l);
    end
end
a_3  = sum(sum(A));
end
function [a_3k]       = Aa3k(lambdaR,lambdaA,epsilon)
global NGroup zeta2 C
for k=1:NGroup
    for l=1:NGroup
        alpha     = C(k,l)*(1/(lambdaA(k,l)-3)-1/(lambdaR(k,l)-3));
        F4        = Falpha(alpha,4);
        F5        = Falpha(alpha,5);
        F6        = Falpha(alpha,6);
        a_3k(k,l) = -epsilon(k,l)^3*F4*zeta2*exp(F5*zeta2+F6*zeta2^2);
        l         = l+1;
    end
    k=k+1;
end
end
%% Association Term
function [A_A]        = Assoc(xComp,GMultiplicity,NParticles,n_Assoc,BondVol,AssocEpsilon,Temp,Vol,lambdaR,epsilon,sigma)
global NComponent NGroup NSites k_B
if sum(sum(sum(AssocEpsilon)))>0
    x0 = zeros(NComponent,NGroup,NSites);
else
    x0 = ones(NComponent,NGroup,NSites);
end
options = optimoptions(@fsolve,'Display','none','OptimalityTolerance',1e-200,'TolFun',1e-30,'TolX',1e-30,'Algorithm','levenberg-marquardt');
F       = @(X_Assoc) AssociationFraction(NParticles,xComp,GMultiplicity,n_Assoc,AssocEpsilon,BondVol,Temp,Vol,X_Assoc,lambdaR,epsilon,sigma);
X_Assoc = fsolve(F,x0,options);
for i=1:NComponent
    for k=1:NGroup
        for a=1:NSites
            A(a) = n_Assoc(a,k)*(log(X_Assoc(i,k,a))+(1-X_Assoc(i,k,a))/2);
            a    = a+1;
        end
        B(k) = GMultiplicity(k,i)*sum(A);
        k    = k+1;
    end
    C(i) = xComp(i)*sum(B);
end
A_A = NParticles*k_B*Temp*sum(C);
end
function [F]          = AssociationFraction(NParticles,xComp,GMultiplicity,n_Assoc,AssocEpsilon,BondVol,Temp,Vol,X_Assoc,lambdaR,epsilon,sigma)
global N_A NComponent NGroup NSites
Delta=AssocStrength(AssocEpsilon,BondVol,NParticles,Vol,Temp,lambdaR,epsilon,sigma);
for i=1:NComponent
    for k=1:NGroup
        for a=1:NSites
            for j=1:NComponent
                for l=1:NGroup
                    for b=1:NSites
                        A(b) = NParticles*xComp(j)*GMultiplicity(l,j)*n_Assoc(b,l)*X_Assoc(j,l,b)*Delta(i,j,k,l,a,b)/Vol;
                        b    = b+1;
                    end
                    B(l) = sum(A);
                    l    = l+1;
                end
                C(j) = sum(B);
                j    = j+1;
            end
            F(i,k,a) = 1/(1+sum(C))-X_Assoc(i,k,a);
            a        = a+1;
        end
        k = k+1;
    end
    i = i+1;
end
end
function [Delta]      = AssocStrength(AssocEpsilon,BondVol,NParticles,Vol,Temp,lambdaR,epsilon,sigma)
global NComponent NGroup NSites zetaF EffD3 EffEpsilon k_B EffSigma xSegment
AssocCorr=[0.075642518,-0.128667137,0.128350632,-0.072532178,0.025778255,-0.006011701,0.000933363,-9.55607e-05,6.19576e-06,-2.30467e-07,3.74606e-09;0.134228218,-0.182682169,0.077166241,-0.000717459,-0.008724273,0.002979718,-0.000484864,4.35262e-05,-2.07789e-06,4.13749e-08,0;-0.565116429,1.009306922,-0.660166946,0.214492212,-0.038846299,0.00406017,-0.000239516,7.25488e-06,-8.58905e-08,0,0;-0.387336383,-0.21161457,0.450442894,-0.176931753,0.031717152,-0.002913689,0.000130194,-2.14506e-06,0,0,0;2.137131809,-2.027984601,0.336709256,0.001181065,-0.006000584,0.000626344,-2.03636e-05,0,0,0,0;-0.300527495,2.899207145,-0.56713484,0.051808513,-0.002393268,4.15107e-05,0,0,0,0,0;-6.210280657,-1.928833603,0.284109761,-0.015760677,0.000368599,0,0,0,0,0,0;11.60835328,0.742215545,-0.082397653,0.001861677,0,0,0,0,0,0,0;-10.26325355,-0.125035689,0.011429914,0,0,0,0,0,0,0,0;4.652974468,-0.001925181,0,0,0,0,0,0,0,0,0;-0.86729622,0,0,0,0,0,0,0,0,0,0];
for k=1:NGroup
    for l=1:NGroup
        A(l) = xSegment(k)*xSegment(l)*sigma(k,l)^3;
    end
    B(k) = sum(A);
end
SigmaX3 = sum(B);
for i=1:NComponent
    for k=1:NGroup
        for a=1:NSites
            for j=1:NComponent
                EffE(i,j) = sqrt(EffSigma(i,i)^3*EffSigma(j,j)^3)*sqrt(EffEpsilon(i,i)*EffEpsilon(j,j))/EffSigma(i,j)^3;
                %gHSd(i,j) = 1/(1-zetaF(4))+3*(EffD3(i,i)*EffD3(j,j))^(1/3)*zetaF(3)/((1-zetaF(4))^2*(EffD3(i,i)^(1/3)+EffD3(j,j)^(1/3)))+2*(EffD3(i,i)*EffD3(j,j))^(2/3)*zetaF(3)^2/((1-zetaF(4))^3*(EffD3(i,i)^(1/3)+EffD3(j,j)^(1/3))^2);
                for l=1:NGroup
                    for n=1:11
                        A = zeros(12-n,1);
                        for m=1:12-n
                            A(m) = AssocCorr(n,m)*(NParticles*SigmaX3/Vol)^(n-1)*(Temp*k_B/EffE(i,j))^(m-1);
                        end
                        B(n) = sum(A);
                    end
                    I(k,l) = sum(B);
                    for b=1:NSites
                        Delta(i,j,k,l,a,b) = (exp(AssocEpsilon(k,l,a,b)/(k_B*Temp))-1)*BondVol(k,l,a,b)*I(k,l);
                        b                  = b+1;
                    end
                    l = l+1;
                end
                j = j+1;
            end
            a = a+1;
        end
        k = k+1;
    end
    i = i+1;
end
end
%% Chain Term
function [A_C]        = Chain(lambdaR,lambdaA,epsilon,sigma,GMultiplicity,NSegment,ShapeFactor,Temp,xComp,NParticles,Vol)
global z EffD3 zetaC NComponent NGroup k_B EffSigma rho_s EffEpsilon

z          = GroupFracComp(GMultiplicity,NSegment,ShapeFactor);
EffD3      = EffectiveDiameter();
EffSigma   = EffectiveSigma(sigma);
EffLambdaR = EffectiveLambda(lambdaR);
EffLambdaA = EffectiveLambda(lambdaA);
EffEpsilon = EffectiveEpsilon(epsilon);
zetaC      = zetaChain(xComp,GMultiplicity,NSegment,ShapeFactor,rho_s);
gMie       = PDFMie(EffSigma,EffLambdaR,EffLambdaA,EffEpsilon,Temp,xComp,GMultiplicity,NSegment,ShapeFactor,Vol);
for i=1:NComponent
    for k=1:NGroup
        if GMultiplicity(k,i)*NSegment(k)*ShapeFactor(k)==1
        A(k) = 0;
        k    = k+1;
        else
        A(k) = (GMultiplicity(k,i)*NSegment(k)*ShapeFactor(k)-1)*log(gMie(i,i));
        k    = k+1;
        end
    end
    B(i) = xComp(i)*sum(A);
    i    = i+1;
end

A_C = -NParticles*k_B*Temp*sum(B);

end
function [EffLambda]  = EffectiveLambda(lambda)
global NComponent NGroup z
for i=1:NComponent
    for j=1:NComponent
        for k=1:NGroup
            for l=1:NGroup
                A(k,l) = z(k,i)*z(l,j)*lambda(k,l);
                l      = l+1;
            end
            k = k+1;
        end
        EffLambda(i,j) = sum(sum(A));
        j              = j+1;
    end
    i = i+1;
end
end
function [EffEpsilon] = EffectiveEpsilon(epsilon)
global NComponent NGroup z
for i=1:NComponent
    for j=1:NComponent
        for k=1:NGroup
            for l=1:NGroup
                A(k,l) = z(k,i)*z(l,j)*epsilon(k,l);
                l      = l+1;
            end
            k = k+1;
        end
        EffEpsilon(i,j) = sum(sum(A));
        j               = j+1;
    end
    i = i+1;
end
end
function [gMie]       = PDFMie(EffSigma,EffLambdaR,EffLambdaA,EffEpsilon,Temp,xComp,GMultiplicity,NSegment,ShapeFactor,Vol)
global k_B NComponent
gHS = PDFHS(EffSigma,Vol);
g1  = PDF1(EffSigma,EffEpsilon,EffLambdaR,EffLambdaA,xComp,GMultiplicity,NSegment,ShapeFactor,Vol);
g2  = PDF2(EffEpsilon,EffLambdaR,EffLambdaA,Temp,xComp,GMultiplicity,NSegment,ShapeFactor,Vol);

for i=1:NComponent
    gMie(i,i) = gHS(i,i)*exp((EffEpsilon(i,i)*g1(i,i)/(k_B*Temp*gHS(i,i)))+(EffEpsilon(i,i)^2*g2(i,i)/((k_B*Temp)^2*gHS(i,i))));
    i         = i+1;
end
end
function [gHS]        = PDFHS(EffSigma,Vol)
global zetaC NComponent EffD3 xEff rho_s

k0 = -log(1-zetaC)+(42*zetaC-39*zetaC^2+9*zetaC^3-2*zetaC^4)/(6*(1-zetaC)^3);
k1 = (zetaC^4+6*zetaC^2-12*zetaC)/(2*(1-zetaC)^3);
k2 = -3*zetaC^2/(8*(1-zetaC)^2);
k3 = (-zetaC^4+3*zetaC^2+3*zetaC)/(6*(1-zetaC)^3);

for i=1:NComponent
    xEff(i,i) = EffSigma(i,i)/(EffD3(i,i))^(1/3);
    gHS(i,i)  = exp(k0+k1*xEff(i,i)+k2*xEff(i,i)^2+k3*xEff(i,i)^3);
    i         = i+1;
end
end
function [g1]         = PDF1(EffSigma,EffEpsilon,EffLambdaR,EffLambdaA,xComp,GMultiplicity,NSegment,ShapeFactor,Vol)
global NComponent EffD3 xEff EffC rho_s

EffC      = EffectiveC(EffLambdaR,EffLambdaA);
Effa1sA   = EffectiveA1s(EffEpsilon,EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Effa1sR   = EffectiveA1s(EffEpsilon,EffLambdaR,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBkA    = EffectiveBk(EffEpsilon,EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBkR    = EffectiveBk(EffEpsilon,EffLambdaR,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);

for j=1:NComponent
    Effa1(j,j) = EffC(j,j)*(xEff(j,j)^EffLambdaA(j,j)*(Effa1sA(j,j)+EffBkA(j,j))-xEff(j,j)^EffLambdaR(j,j)*(Effa1sR(j,j)+EffBkR(j,j)));
    j          = j+1;
end
Da1Drho_s = PartialA1(EffLambdaA,EffLambdaR,EffEpsilon);

for i=1:NComponent
    xEff(i,i) = EffSigma(i,i)/(EffD3(i,i))^(1/3);
    g1(i,i)   = (1/(2*pi*EffEpsilon(i,i)*EffD3(i,i)))*(3*Da1Drho_s(i,i)-EffC(i,i)*EffLambdaA(i,i)*xEff(i,i)^EffLambdaA(i,i)*(Effa1sA(i,i)+EffBkA(i,i))/rho_s+EffC(i,i)*EffLambdaR(i,i)*xEff(i,i)^EffLambdaR(i,i)*(Effa1sR(i,i)+EffBkR(i,i))/rho_s);
    i         = i+1;
end
end
function [Effa1s]     = EffectiveA1s(EffEpsilon,EffLambda,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor)
global NComponent EffD3 
zetaC      = zetaChain(xComp,GMultiplicity,NSegment,ShapeFactor,rho_s);
for i=1:NComponent
    A             = [ 0.81096   1.7888  -37.578   92.284;
                      1.02050 -19.3410  151.260 -463.500; 
                     -1.90570  22.8450 -228.140  973.920; 
                      1.08850  -6.1962  106.980 -677.640 ]*[1; 1/EffLambda(i,i); 1/EffLambda(i,i)^2; 1/EffLambda(i,i)^3];
    zeta_eff(i,i) =   A(1).*zetaC+A(2).*zetaC^2+A(3).*zetaC^3+A(4).*zetaC^4;
    Effa1s(i,i)   =  -2*pi*rho_s*EffEpsilon(i,i)*EffD3(i,i)*(1-zeta_eff(i,i)/2)/((EffLambda(i,i)-3)*(1-zeta_eff(i,i))^3);
    i             =   i+1;
end
end
function [EffBk]      = EffectiveBk(EffEpsilon,EffLambda,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor)
global xEff NComponent EffD3
zetaC      = zetaChain(xComp,GMultiplicity,NSegment,ShapeFactor,rho_s);
for i=1:NComponent
    I          = (1-xEff(i,i).^(3-EffLambda(i,i)))/(EffLambda(i,i)-3);
    J          = (1-xEff(i,i).^(4-EffLambda(i,i))*(EffLambda(i,i)-3)+xEff(i,i)^(3-EffLambda(i,i))*(EffLambda(i,i)-4))/((EffLambda(i,i)-4)*(EffLambda(i,i)-3));
    EffBk(i,i) = 2*pi*rho_s*EffD3(i,i)*EffEpsilon(i,i)*((1-zetaC/2)*I/(1-zetaC)^3-(9*zetaC*(1+zetaC)*J)/(2*(1-zetaC)^3));
    i          = i+1;
end
end
function [Da1Drho_s]  = PartialA1(EffLambdaA,EffLambdaR,EffEpsilon)
global NComponent EffC xEff rho_s 
Da1sDrho_sA = PartialA1s(EffLambdaA,EffEpsilon);
DBkDrho_sA  = PartialBk(EffLambdaA,EffEpsilon);
Da1sDrho_sR = PartialA1s(EffLambdaR,EffEpsilon);
DBkDrho_sR  = PartialBk(EffLambdaR,EffEpsilon);
for i=1:NComponent
    Da1Drho_s(i,i) = EffC(i,i)*(xEff(i,i)^EffLambdaA(i,i)*(Da1sDrho_sA(i,i)+DBkDrho_sA(i,i))-xEff(i,i)^EffLambdaR(i,i)*(Da1sDrho_sR(i,i)+DBkDrho_sR(i,i)));
    i              = i+1;
end
end
function [Da1sDrho_s] = PartialA1s(EffLambda,EffEpsilon)
global NComponent zeta1 rho_s EffD3 
for i=1:NComponent
    A               = [ 0.81096   1.7888  -37.578   92.284;
                        1.02050 -19.3410  151.260 -463.500; 
                       -1.90570  22.8450 -228.140  973.920; 
                        1.08850  -6.1962  106.980 -677.640 ]*[1; 1/EffLambda(i,i); 1/EffLambda(i,i)^2; 1/EffLambda(i,i)^3];
    zetaEff(i,i)    = A(1).*zeta1+A(2).*zeta1^2+A(3).*zeta1^3+A(4).*zeta1^4;
    DzetaEff(i,i)   = (A(1)+2*A(2).*zeta1+3*A(3).*zeta1^2+4*A(4).*zeta1^3).*zeta1/rho_s;
    Da1sDrho_s(i,i) = -2*pi*(EffEpsilon(i,i)*EffD3(i,i))/(EffLambda(i,i)-3)*((1-zetaEff(i,i)/2)/(1-zetaEff(i,i))^3+rho_s*((3*(1-zetaEff(i,i)/2)*(1-zetaEff(i,i))^2-0.5*(1-zetaEff(i,i))^3)/(1-zetaEff(i,i))^6)*DzetaEff(i,i));
    i               = i+1;
end
end
function [DBkDrho_s]  = PartialBk(EffLambda,EffEpsilon)
global NComponent EffD3  xEff zeta1 rho_s

for i=1:NComponent
    I              = (1-xEff(i,i).^(3-EffLambda(i,i)))/(EffLambda(i,i)-3);
    J              = (1-xEff(i,i).^(4-EffLambda(i,i))*(EffLambda(i,i)-3)+xEff(i,i)^(3-EffLambda(i,i))*(EffLambda(i,i)-4))/((EffLambda(i,i)-4)*(EffLambda(i,i)-3));
    DBkDrho_s(i,i) = 2*pi*EffD3(i,i)*EffEpsilon(i,i)*(((1-zeta1/2)*I/(1-zeta1)^3-9*zeta1*(1+zeta1)*J/(2*(1-zeta1)^3))+rho_s*((3*(1-zeta1/2)*(1-zeta1)^2-0.5*(1-zeta1)^3)*I/(1-zeta1)^6-9*J*((1+2*zeta1)*(1-zeta1)^3+zeta1*(1+zeta1)*3*(1-zeta1)^2)/(2*(1-zeta1)^6))*zeta1/rho_s);
    i              = i+1;
end
end
function [g2]         = PDF2(EffEpsilon,EffLambdaR,EffLambdaA,Temp,xComp,GMultiplicity,NSegment,ShapeFactor,Vol)
global NComponent rho_s

yC   = EffFac(EffEpsilon,EffLambdaR,EffLambdaA,Temp);
gMCA = PDFMCA(EffEpsilon,EffLambdaR,EffLambdaA,xComp,GMultiplicity,NSegment,ShapeFactor);

for i=1:NComponent
    g2(i,i) = (1+yC(i,i))*gMCA(i,i);
    i       = i+1;
end
end
function [yC]         = EffFac(EffEpsilon,EffLambdaR,EffLambdaA,Temp)
global NComponent zeta2 k_B EffC EffAlpha
Phi=[7.5365557,-359.4400,1550.9,-1.19932,-1911.28,9236.9,10;-37.60463,1825.6,-5070.1,9.063632,21390.175,-129430,10;71.745953,-3168,6534.6,-17.9482,-51320.7,357230,0.57;-46.83552,1884.2,-3288.7,11.34027,37064.54,-315530,-6.7;-2.467982,-0.82376,-2.7171,20.52142,1103.742,1390.2,-8;-0.50272,-3.1935,2.0883,-56.6377,-3264.61,-4518.2,NaN;8.0956883,3.7090,0,40.53683,2556.181,4241.6,NaN];
for i=1:NComponent
    EffAlpha(i,i) = EffC(i,i)*(1/(EffLambdaA(i,i)-3)-1/(EffLambdaR(i,i)-3));
    yC(i,i)       = Phi(1,7)*(-tanh(Phi(2,7)*(Phi(3,7)-EffAlpha(i,i)))+1)*zeta2*(exp(EffEpsilon(i,i)/(k_B*Temp))-1)*exp(Phi(4,7)*zeta2+Phi(5,7)*zeta2^2);
    i             = i+1;
end
end
function [gMCA]       = PDFMCA(EffEpsilon,EffLambdaR,EffLambdaA,xComp,GMultiplicity,NSegment,ShapeFactor)
global NComponent K_HS EffD3 EffC xEff rho_s
Effa1s2A  = EffectiveA1s(EffEpsilon,2*EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Effa1s2R  = EffectiveA1s(EffEpsilon,2*EffLambdaR,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBk2A   = EffectiveBk(EffEpsilon,2*EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBk2R   = EffectiveBk(EffEpsilon,2*EffLambdaR,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Effa1sAR  = EffectiveA1s(EffEpsilon,EffLambdaR+EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBkAR   = EffectiveBk(EffEpsilon,EffLambdaR+EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Da2Drho_s = PartialA2(EffLambdaA,EffLambdaR,EffEpsilon,xComp,GMultiplicity,NSegment,ShapeFactor);

for i=1:NComponent
    gMCA(i,i) = 1/(2*pi*EffEpsilon(i,i)^2*EffD3(i,i))*(3*Da2Drho_s(i,i)-EffEpsilon(i,i)*K_HS*EffC(i,i)^2*EffLambdaR(i,i)*xEff(i,i)^(2*EffLambdaR(i,i))*(Effa1s2R(i,i)+EffBk2R(i,i))/rho_s+EffEpsilon(i,i)*K_HS*EffC(i,i)^2*(EffLambdaR(i,i)+EffLambdaA(i,i))*xEff(i,i)^(EffLambdaR(i,i)+EffLambdaA(i,i))*(Effa1sAR(i,i)+EffBkAR(i,i))/rho_s-EffEpsilon(i,i)*K_HS*EffC(i,i)^2*EffLambdaA(i,i)*xEff(i,i)^(2*EffLambdaA(i,i))*(Effa1s2A(i,i)+EffBk2A(i,i))/rho_s);
    i         = i+1;
end
end
function [Da2Drho_s]  = PartialA2(EffLambdaA,EffLambdaR,EffEpsilon,xComp,GMultiplicity,NSegment,ShapeFactor)
global NComponent EffC zeta1 xEff rho_s
DK_HS        = -((4*(1-zeta1)^3*(1+4*zeta1+4*zeta1^2-4*zeta1^3+zeta1^4)+(1-zeta1)^4*(4+8*zeta1-12*zeta1^2+4*zeta1^3))/(1+4*zeta1+4*zeta1^2-4*zeta1^3+zeta1^4)^2)*zeta1/rho_s;
K_HS         = (1-zeta1)^4/(1+4*zeta1+4*zeta1^2-4*zeta1^3+zeta1^4);
Effa1s2A     = EffectiveA1s(EffEpsilon,2*EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Effa1s2R     = EffectiveA1s(EffEpsilon,2*EffLambdaR,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBk2A      = EffectiveBk(EffEpsilon,2*EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBk2R      = EffectiveBk(EffEpsilon,2*EffLambdaR,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Effa1sAR     = EffectiveA1s(EffEpsilon,EffLambdaR+EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
EffBkAR      = EffectiveBk(EffEpsilon,EffLambdaR+EffLambdaA,rho_s,xComp,GMultiplicity,NSegment,ShapeFactor);
Da1s2ADrho_s = PartialA1s(2*EffLambdaA,EffEpsilon);
Da1s2RDrho_s = PartialA1s(2*EffLambdaR,EffEpsilon);
Da1sARDrho_s = PartialA1s(EffLambdaA+EffLambdaR,EffEpsilon);
DBk2ADrho_s  = PartialBk(2*EffLambdaA,EffEpsilon);
DBk2RDrho_s  = PartialBk(2*EffLambdaR,EffEpsilon);
DBkARDrho_s  = PartialBk(EffLambdaA+EffLambdaR,EffEpsilon);

for i=1:NComponent
    Da2Drho_s(i,i) = 0.5*EffEpsilon(i,i)*EffC(i,i)^2*(DK_HS*(xEff(i,i)^(2*EffLambdaA(i,i))*(Effa1s2A(i,i)+EffBk2A(i,i))-2*xEff(i,i)^(EffLambdaA(i,i)+EffLambdaR(i,i))*(Effa1sAR(i,i)+EffBkAR(i,i))+xEff(i,i)^(2*EffLambdaR(i,i))*(Effa1s2R(i,i)+EffBk2R(i,i)))+K_HS*(xEff(i,i)^(2*EffLambdaA(i,i))*(Da1s2ADrho_s(i,i)+DBk2ADrho_s(i,i))-2*xEff(i,i)^(EffLambdaA(i,i)+EffLambdaR(i,i))*(Da1sARDrho_s(i,i)+DBkARDrho_s(i,i))+xEff(i,i)^(2*EffLambdaR(i,i))*(Da1s2RDrho_s(i,i)+DBk2RDrho_s(i,i))));
    i              = i+1;
end
end