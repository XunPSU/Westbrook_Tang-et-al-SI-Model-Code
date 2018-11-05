function dx= STARModel(t,x,par)

dx = zeros(5,1);
%%%% this model is the same as the IFFL model


S=x(1);  %%% STAR 
Y=x(2);  %%% mRNA
Yi=x(3);  
G=x(4);
Gm=x(5);  %%% measurable GFP

Py=par.Pytot;  %%% free GFP plasmid

dx(1) = par.alpha_s*par.Ps - par.deg_s*S - par.beta_s*S*Py;  %%% STAR

dx(2) = par.alpha_m*par.beta_s*S*Py - par.deg_m*Y - par.KI*Y + par.KE*Yi; %%% mRFP production

dx(3) = par.KI*Y - par.KE*Yi; %%% maturation

dx(4) = par.KE*Yi - par.alpha_gm*G; %%% elongation 

dx(5) = par.alpha_gm*G; %%%measurement

return 