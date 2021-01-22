%% derive analytical solution for Cauchy stress for hydrostatic testing

%%
clear all;
clc;

%%% define symbolic parameters
syms lam mu;
syms A B C;
syms e1;
assume(e1,'real');

%%% functions for the second and third-order elastic constants
I = @(i,j,k,l) ((i==k)*(j==l)+(i==l)*(j==k))/2;  % Eq.(6)

SC = @(i,j,k,l) lam*(i==j)*(k==l) + 2*mu*I(i,j,k,l); % Eq.(4)

TC = @(i,j,k,l,m,n) 2*C*(i==j)*(k==l)*(m==n) + ...
    2*B*((i==j)*I(k,l,m,n) + (k==l)*I(m,n,i,j) + (m==n)*I(i,j,k,l)) + ...
    A/2*((i==k)*I(j,l,m,n) + (i==l)*I(j,k,m,n) + (j==k)*I(i,l,m,n) + (j==l)*I(i,k,m,n));  %Eq.(5)


%%% deformation matrix & Green-Lagrangian strain tensor
F = [e1+1 0 0; 0 e1+1 0; 0 0 e1+1]; %Eq. (8)
EE = 0.5*(F'*F-eye(3)); % Eq.(9)
J=det(F); % determinant of F

%%% Loop over repeated indices Eq.(3)
S = sym(zeros(3,3)); % prellocate location to store PK-II stress
for i= 1:3
    for j=1:3

        part1 = 0; % second-order part
        part2 = 0; % third-order part

        for k= 1:3
            for l = 1:3
                part1 = part1 + SC(i,j,k,l)*EE(k,l);
                
                for m = 1:3
                    for n = 1:3
                        part2 = part2 + 1/2*TC(i,j,k,l,m,n)*EE(k,l)*EE(m,n);
                    end
                end
                           
            end
        end
        
       S(i,j)= expand(part1) + expand(part2); % store into S 
               
    end
end

%%% change PK-II stress S to the cauchy stress
Cauchy = F*S*F'/J ; 

%%% Taylor expansion of stress to the quadratic of e1
%Stress11 = taylor(Cauchy(1,1),e1,0,'Order',3);

Stress11 = collect(Cauchy(1,1),[A,B,C,lam,mu]);
pretty(Stress11) % then ignore higher order term (third or higher) then get Eq.(10)

% %% validate with the numerical outputs
% lam = 4.45*10^9;
% mu = 4.13*10^9;
% A = -5011*10^9;
% B = -1692*10^9;
% C = -158*10^9;
% e1 = (-0.0006:0.00001:0);
% cal_sigxx = (e1.^2*A+9*e1.^2*B+9*e1.^2*C+(1.5*e1.^2+3*e1)*lam+mu*(e1.^2+2*e1))./(e1+1); % Eq. (10)
% plot(-e1,-cal_sigxx); % in the equation, we assume the compression is negative. 
% hold on; % but in t he plot, we add negation to better illustrate. 
% 
% % input output from comsol
% file = textread('outfromcomsol_hydr_noheader.txt');
% ee1 = file(:,2);
% ss1 = file(:,4);
% plot(-ee1,-ss1,'*');
% hold off;











