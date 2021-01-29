%% derive analytical solution for PK-I stress based on the linear strain e1,e2,e3 for general loading

%%
clear all;
clc;

%%% define symbolic parameters
syms lam mu; 
syms A B C;
syms e1 e2 e3;
assume(e1,'real');
assume(e2,'real');
assume(e3,'real');

%%% functions for the second and third-order elastic constants
I = @(i,j,k,l) ((i==k)*(j==l)+(i==l)*(j==k))/2;  % fourth order identity

SC = @(i,j,k,l) lam*(i==j)*(k==l) + 2*mu*I(i,j,k,l); % Eq.(4)

TC = @(i,j,k,l,m,n) 2*C*(i==j)*(k==l)*(m==n) + ...
    2*B*((i==j)*I(k,l,m,n) + (k==l)*I(m,n,i,j) + (m==n)*I(i,j,k,l)) + ...
    A/2*((i==k)*I(j,l,m,n) + (i==l)*I(j,k,m,n) + (j==k)*I(i,l,m,n) + (j==l)*I(i,k,m,n));  %Eq.(5)


%%% deformation matrix & Green-Lagrangian strain tensor
F = [e1+1 0 0; 0 e2+1 0; 0 0 e3+1]; %Eq. (7)
EE = 0.5*(F'*F-eye(3)); % Eq.(8)
%J=det(F); % determinant of F

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

%%% change PK-II stress S to the PK-I stress P
PK_I = expand(F*S) ;  % here we get the stress tensor. we only cares about normal stresses. No shear stresses. 

%%
%%% Collect based on A,B,C,lam,mu
%collA = coeffs(Cauchy(1,1),A); % collect based on A using polonominals 
%coffA = collA(1,2); % the coefficient of A 
%fprintf('For A coefficient:\n %s',char(coffA))

P11 = collect(PK_I(1,1),[A,B,C,lam,mu]);
pretty(P11) % ignore third- or higher terms and get Eq. (9)

P22 = collect(PK_I(2,2),[A,B,C,lam,mu]);
pretty(P22) % ignore third- or higher terms and get Eq. (10)

P33 = collect(PK_I(3,3),[A,B,C,lam,mu]);
pretty(P33) % ignore third- or higher terms and get Eq. (10)




