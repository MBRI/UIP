%============================================================%
% Solve a cubic equation in SolvingCubicEquation.pdf         %
% Date: January 6, 2015                                      %                
%============================================================%

% Risk premium depends on nominal interest rate
% Lagged interest rate in Taylor rule
% Shock is to Phillilps curve..equilibrium real interest rate
% Does allow serial correlation in risk premium shock
% Two shocks are uncorrelated
% Output written to SolvingCubicEquationsNagelPhillips.xls
% Notes from 1/6/15

%---------------------%
%Try different values
%---------------------%
%known parameter values 
bet = 0.998;
rho = 0.915;
del = 0.0141;
sig = 1.5*(1-rho);
alp = .15;
psi = 0.99;
mu = 0;
var_qbar = 1;
var_eta = 0.04;
;

syms x0
%Solve the following equation to find the values of c
c0 = double(solve(del*x0*(1+alp)*sig-del*(x0-rho)*x0-(x0-1)*(1-bet*x0)*(x0-rho)==0));

%Keep the one that is less than one in absolute value
c = c0(find(lt(abs(c0),1)==1));

c1 = c0(find(gt(abs(c0),1)==1));

%With the solution for c, solve for other variables
k = (c-rho)/sig;
f = k*(1-bet*c)/del;
a = (psi-1)*sig*del/((1-psi)*(1-bet*(psi+k*sig))+del*(sig*(1+alp-f-k)-psi));
d = a*(sig*(1+alp-f-k)-psi)/(sig*(psi-1));
g = a/sig;
b = sig*del/((1-mu)*(1-bet*mu-sig*bet*k)-del*mu-sig*del*(f+k-1-alp));
h = b/sig;
e = (1-bet*mu-sig*bet*k)*b/(sig*del);
m = a*(1-k)-g*psi;
n = b*(1-k)-h*mu;
p = c*(1-k);
q = alp*a;
r = alp*b-1;
s = alp*c;
v = q/((1-c)*(1-psi));
w = (1/(1-mu))*((r+1)/(1-c)-1);
x = s/(1-c);
covqbarnom=(a*psi/(1-psi*c))*var_qbar;
covetanom=(b*mu/(1-mu*c))*var_eta;
varnom = ((a^2)/(1-c^2))*var_qbar+((b^2)/(1-c^2))*var_eta+(2*a*c/(1-c^2))*covqbarnom+(2*b*c/(1-c^2))*covetanom;
corlam = m*q*var_qbar+n*r*var_eta+p*s*varnom+(m*s+q*p)*covqbarnom+(n*s+r*p)*covetanom;
corblam = m*v*var_qbar+n*w*var_eta+p*x*varnom+(m*x+v*p)*covqbarnom+(n*x+w*p)*covetanom;
varreal = (m^2)*var_qbar+(n^2)*var_eta+(p^2)*varnom+2*m*p*covqbarnom+2*n*p*covetanom;
varlam = (q^2)*var_qbar+(r^2)*var_eta+(s^2)*varnom+2*q*s*covqbarnom+2*r*s*covetanom;
betshort = corlam/varreal;
betlong = corblam/varreal;
varinf = (g^2)*var_qbar+(h^2)*var_eta+(k^2)*varnom+2*g*k*covqbarnom+2*h*k*covetanom;
sercov = g*(g*psi+k*a)*var_qbar+h*(h*mu+k*b)*var_eta+(k^2)*c*varnom+k*(k*a+g*(psi+c))*covqbarnom+k*(k*b+h*(mu+c))*covetanom;
sercorinf=sercov/varinf;
sernom = a*covqbarnom+b*covetanom+c*varnom;
sercorint = sernom/varnom;
covrealnom = a*m*var_qbar+b*n*var_eta+c*p*varnom+(a*p+c*m)*covqbarnom+(b*p+c*n)*covetanom;

% Wrap up parameters and variables. 
d1 = [{'sigma', 'rho', 'alpha', 'beta', 'delta', 'var_qbar', 'var_eta', 'psi', 'mu', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'k', 'm', 'n', 'p', 'q', 'r', 's', 'v', 'w', 'x', 'corlam', 'corblam', 'betshort', 'betlong', 'varreal', 'varlam', 'varnom', 'varinf' 'sercorinf' 'sercorint' 'covrealnom' }; 
    [num2cell(sig) num2cell(rho) num2cell(alp) num2cell(bet) num2cell(del) num2cell(var_qbar) num2cell(var_eta) num2cell(psi) num2cell(mu) num2cell(a) num2cell(b) num2cell(c) num2cell(d) num2cell(e) num2cell(f) num2cell(g) num2cell(h) num2cell(k) num2cell(m) num2cell(n) num2cell(p) num2cell(q) num2cell(r) num2cell(s) num2cell(v) num2cell(w) num2cell(x) num2cell(corlam) num2cell(corblam) num2cell(betshort) num2cell(betlong) num2cell(varreal) num2cell(varlam) num2cell(varnom) num2cell(varinf) num2cell(sercorinf) num2cell(sercorint) num2cell(covrealnom) ]]'; 

%Write the results into a Excel Spreadsheet
xlswrite('SolvingCubicEquationsNagelPhillipsPersistence.xls', d1, 'results', 'A1');
