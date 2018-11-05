%%%% Pulse model fitting
close all
clear all


Con2=4;

%%%% Experimental results of Pulse network

M10=csvread('Data3010.csv',14,0);  %%%% load Pulse experimental measurements
Exp=M10(:,29:37);
Cexp=mean(Exp,2);

%%%% paramter constraints

Lalpha=10^(-2);
Ualpha=10;
stepRalpha=(Ualpha-Lalpha)/1000; 

LaM=10^(-1); 
UaM=100;
stepRaM=(UaM-LaM)/1000; 

Lgamma=10^3;
Ugamma=10^7;
stepGamma=(Ugamma-Lgamma)/1000; 

Ldeg=10^(-5);
Udeg=10^(-1);
stepDeg=(Udeg-Ldeg)/1000; 

Lomega=10^3;
Uomega=10^7;
stepOmega=(Uomega-Lomega)/1000;  

Lbd=10^3;
Ubd=10^7;
stepBeta=(Ubd-Lbd)/1000; 

LKi=10^(-4);
UKi=10^(-2);
stepKi=(UKi-LKi)/1000; 

LKe=10^(-4);
UKe=10^(-2);
stepKe=(UKe-LKe)/1000; 

LGm=10^(-3);
UGm=10^(-1);
stepGm=(UGm-LGm)/1000; 


Percent=[0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95];

for nj = 1:10
figure(nj)    
    savefile=strcat('IFFLFit105000Final',int2str(nj),'.mat');
    Sample=randi([1 10],1,15);  
    
    %%%% crispr parameters
    par.alpha_cr=Lalpha+(Ualpha-Lalpha)*Percent(Sample(1));
    par.deg_cr=Ldeg+(Udeg-Ldeg)*Percent(Sample(2));
    par.deg_tr=Ldeg+(Udeg-Ldeg)*Percent(Sample(3));
    par.deg_g=Ldeg+(Udeg-Ldeg)*Percent(Sample(4));
    
    par.omega=Lomega+(Uomega-Lomega)*Percent(Sample(5));
    par.gamma1=Lgamma+(Ugamma-Lgamma)*Percent(Sample(6));  
    par.gamma2=Lgamma+(Ugamma-Lgamma)*Percent(Sample(7)); 
    
    %%% star paramters
    par.alpha_s=Lalpha+(Ualpha-Lalpha)*Percent(Sample(8));
    par.beta_s=Lbd+(Ubd-Lbd)*Percent(Sample(9));
    par.deg_s=Ldeg+(Udeg-Ldeg)*Percent(Sample(10));
     
    %%% EGFP parameters
    par.alpha_m=LaM+(UaM-LaM)*Percent(Sample(11));
    par.deg_m=Ldeg+(Udeg-Ldeg)*Percent(Sample(12));
    par.KI=LKi+(UKi-LKi)*Percent(Sample(13)); 
    par.KE=LKe+(UKe-LKe)*Percent(Sample(14)); 
    par.alpha_gm=LGm+(UGm-LGm)*Percent(Sample(15)); 
   
    
    %%%% exp. settings
    par.Pytot=0.5*10^(-9);
    par.Pcr=0.25*10^(-9);
    par.Ptr=0.25*10^(-9);
    par.Ps=16*10^(-9);      
    par.dCas9tot=35*10^(-9);
    
    tspan=0:300:14400;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    y0=[0 0 0 0 par.dCas9tot 0 0 0 0 0 0];
    [t,y]=ode23s(@PulseModel,tspan,y0,options,par);
    y=y*10^6;
    OError=sum((y(:,11)-Cexp).^2);
    
    iter=1;
    Maxiter=210000;
    en=0;
    Pre = zeros(49,Maxiter);
    SimuResults=zeros(Maxiter,16); 


        while iter<Maxiter  
        clear y Newpar
        Newpar=par;
        step=randi([-1 1],1,15);  
        
        %%% crispr
        Newpar.alpha_cr=par.alpha_cr+step(1)*stepRalpha;
        Newpar.deg_tr=par.deg_tr+step(2)*stepDeg;
        Newpar.deg_cr=par.deg_cr+step(3)*stepDeg;
        Newpar.deg_g=par.deg_g+step(4)*stepDeg;
        
        Newpar.gamma1=par.gamma1+step(5)*stepGamma;
        Newpar.gamma2=par.gamma2+step(6)*stepGamma;
        Newpar.omega=par.omega+step(7)*stepOmega;  
        %%%% star 
        Newpar.alpha_s=par.alpha_s+step(8)*stepRalpha;
        Newpar.deg_s=par.deg_s+step(9)*stepDeg;
        Newpar.beta_s=par.beta_s+step(10)*stepBeta; 
        %%%% EGFP
        Newpar.alpha_m=par.alpha_m+step(11)*stepRaM;
        Newpar.deg_m=par.deg_m+step(12)*stepDeg;
        Newpar.KI=par.KI+step(13)*stepKi; 
        Newpar.KE=par.KE+step(14)*stepKe; 
        Newpar.alpha_gm=par.alpha_gm+step(15)*stepGm; 



        if (Newpar.alpha_cr>=Lalpha) && (Newpar.alpha_cr<=Ualpha)...
            && (Newpar.gamma1>=Lgamma) && (Newpar.gamma1<=Ugamma)...
            && (Newpar.gamma2>=Lgamma) && (Newpar.gamma2<=Ugamma)...
            && (Newpar.deg_tr>=Ldeg) && (Newpar.deg_tr<=Udeg)...
            && (Newpar.deg_cr>=Ldeg) && (Newpar.deg_cr<=Udeg)...
            && (Newpar.deg_g>=Ldeg) && (Newpar.deg_g<=Udeg)...
            && (Newpar.omega>=Lomega) && (Newpar.omega<=Uomega)... 
            && (Newpar.alpha_s>=Lalpha) && (Newpar.alpha_s<=Ualpha)...
            && (Newpar.deg_s>=Ldeg) && (Newpar.deg_s<=Udeg)...
            && (Newpar.beta_s>=Lbd) && (Newpar.beta_s<=Ubd)...
            && (Newpar.alpha_m>=LaM) && (Newpar.alpha_m<=UaM)...
            && (Newpar.deg_m>=Ldeg) && (Newpar.deg_m<=Udeg)...
            && (Newpar.KI>=LKi) && (Newpar.KI<=UKi)...
            && (Newpar.KE>=LKe) && (Newpar.KE<=UKe)...
            && (Newpar.alpha_gm>=LGm) && (Newpar.alpha_gm<=UGm)
        
            en=en+1; 
            
            tspan=0:300:14400; 
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            y0=[0 0 0 0 par.dCas9tot 0 0 0 0 0 0];
            [t,y]=ode23s(@PulseModel,tspan,y0,options,par);
            y=y*10^6;
            Pre(:,iter) = y(:,11);
            NError=sum((y(:,11)-Cexp).^2);
            Pro=exp(-(NError-OError)/(2*0.25^2));  %%% Metropolis criteria with T=0.125 
                                               %%% with standard deviation as 0.25
            if NError<OError
                OError=NError;
                par=Newpar;
                iter=iter+1;
                SimuResults(iter,1)=OError;
            elseif Pro>rand(1)
                OError=NError;
                par=Newpar;
                iter=iter+1;
                SimuResults(iter,1)=OError;
            else
                iter=iter+1;
                SimuResults(iter,1)=OError;
            end
        
        else
        iter=iter+1;
        SimuResults(iter,1)=OError;
        end
        
        SimuResults(iter,2)=par.alpha_s;
        SimuResults(iter,3)=par.deg_s;
        SimuResults(iter,4)=par.beta_s;

        SimuResults(iter,5) = par.alpha_cr;
        SimuResults(iter,6)=par.omega;
        SimuResults(iter,7)=par.gamma1;
        SimuResults(iter,8)=par.gamma2;
        SimuResults(iter,9)=par.deg_tr;
        SimuResults(iter,10)=par.deg_cr;
        SimuResults(iter,11)=par.deg_g;
        
        SimuResults(iter,12)=par.alpha_m;
        SimuResults(iter,13)=par.deg_m;
        SimuResults(iter,14)=par.KI;
        SimuResults(iter,15)=par.KE;
        SimuResults(iter,16)=par.alpha_gm;

        clear y  
        end
 
        
save(savefile,'SimuResults')

[v,loca]=min(SimuResults(2:end,1));
plot(tspan./60,Pre(:,loca+1),'-r',tspan./60,Cexp,'-.r','LineWidth',2)
xlabel('Time (min)')
ylabel('EGFP (\muM)')
legend('Simu','Exp')
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

end

