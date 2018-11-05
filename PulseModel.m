%%%%% Pulse model

function dx= PulseModel(t,x,par)

dx = zeros(11,1);

S=x(1);   
crRNA=x(2);
trRNA=x(3);
gRNA=x(4); 
dCas9=x(5); 
Complex=x(6); 
Pyrep=x(7);  
Y=x(8); 
yi=x(9);  
G=x(10);
Gm=x(11);  

Py=par.Pytot-Pyrep;  


dx(1) = par.alpha_s*par.Ps - par.deg_s*S - par.beta_s*S*Py;  

dx(2) = par.alpha_cr*par.Pcr - par.deg_cr*crRNA - par.gamma1*crRNA*trRNA; 

dx(3) = par.alpha_cr*par.Ptr - par.deg_tr*trRNA - par.gamma1*crRNA*trRNA; 

dx(4) = par.gamma1*crRNA*trRNA - par.gamma2*gRNA*dCas9 - par.deg_g*gRNA; 

dx(5) = - par.gamma2*gRNA*dCas9; 

dx(6) = par.gamma2*gRNA*dCas9 - par.omega*Complex*Py; 

dx(7) = par.omega*Complex*Py; 

dx(8) = par.alpha_m*par.beta_s*S*Py - par.deg_m*Y - par.KI*Y + par.KE*yi; 

dx(9) = par.KI*Y - par.KE*yi; 

dx(10) = par.KE*yi - par.alpha_gm*G;  

dx(11) = par.alpha_gm*G; 

return 