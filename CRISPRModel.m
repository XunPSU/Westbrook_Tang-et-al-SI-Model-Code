%%%%% CRISPR Model
%%%%%
function dx= CRISPRModel(t,x,par)

dx = zeros(10,1);

crRNA=x(1); 
trRNA=x(2);
gRNA=x(3);
dCas9=x(4);  
Complex=x(5); 
Pyrep=x(6); 
Y=x(7); 
yi=x(8); 
G=x(9);
Gm=x(10);  

Py=par.Pytot-Pyrep; 

dx(1) = par.alpha_cr*par.Pcr - par.deg_cr*crRNA - par.gamma1*crRNA*trRNA; 

dx(2) = par.alpha_cr*par.Ptr - par.deg_tr*trRNA - par.gamma1*crRNA*trRNA; 

dx(3) = par.gamma1*crRNA*trRNA - par.gamma2*gRNA*dCas9 - par.deg_g*gRNA; 

dx(4) = - par.gamma2*gRNA*dCas9; 

dx(5) = par.gamma2*gRNA*dCas9 - par.omega*Complex*Py; 

dx(6) = par.omega*Complex*Py; 

dx(7) = par.alpha_m*Py - par.deg_m*Y-par.KI*Y + par.KE*yi; 

dx(8) = par.KI*Y - par.KE*yi; 

dx(9) = par.KE*yi - par.alpha_gm*G; 

dx(10) = par.alpha_gm*G; 


return 