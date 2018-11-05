%%%% Main script for CRISPR model fitting
close all
clear all

Con2=3;
M09=csvread('Data3009.csv',14,0);  %%%% load experimental measurement


%%%% process and preview experimental data
figure
Ave09=zeros(49,Con2);
for i=1:Con2
    for j=1:9
        Ave09(:,i)=Ave09(:,i)+M09(:,(i-1)*9+j+1);
    end
end
Ave09=Ave09./9;


%%%%
%%%%
%%%% paramter constraints

Lalpha=10^(-2);
Ualpha=10;
stepRalpha=(Ualpha-Lalpha)/1000; 

LaM=10^(-1);  
UaM=100;
stepRaM=(UaM-LaM)/1000; 

LKi=10^(-4);
UKi=10^(-2);
stepKi=(UKi-LKi)/1000; 

LKe=10^(-4);
UKe=10^(-2);
stepKe=(UKe-LKe)/1000; 

LGm=10^(-3);
UGm=10^(-1);
stepGm=(UGm-LGm)/1000; 

Lgamma=10^3;
Ugamma=10^7;
stepGamma=(Ugamma-Lgamma)/1000;  

Ldeg=10^(-5);
Udeg=10^(-1);
stepDeg=(Udeg-Ldeg)/1000; 

Lomega=10^3;
Uomega=10^7;
stepOmega=(Uomega-Lomega)/1000;  


Percent=[0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95];

for nj=6:10
    savefile=strcat('M09Fit21000Final',int2str(nj),'.mat');
    
Sample=randi([1 10],1,12);       
par.alpha_cr=Lalpha+(Ualpha-Lalpha)*Percent(Sample(1));
par.alpha_m=LaM+(UaM-LaM)*Percent(Sample(2));
par.deg_cr=Ldeg+(Udeg-Ldeg)*Percent(Sample(3));
par.deg_tr=Ldeg+(Udeg-Ldeg)*Percent(Sample(4));
par.deg_g=Ldeg+(Udeg-Ldeg)*Percent(Sample(5));
par.deg_m=Ldeg+(Udeg-Ldeg)*Percent(Sample(6));
par.omega=Lomega+(Uomega-Lomega)*Percent(Sample(7));
par.gamma1=Lgamma+(Ugamma-Lgamma)*Percent(Sample(8));  
par.gamma2=Lgamma+(Ugamma-Lgamma)*Percent(Sample(9));
par.KI=LKi+(UKi-LKi)*Percent(Sample(10)); 
par.KE=LKe+(UKe-LKe)*Percent(Sample(11)); 
par.alpha_gm=LGm+(UGm-LGm)*Percent(Sample(12));
   

par.Pytot=0.5*10^(-9);  
par.dCas9tot=35*10^(-9);


Simu=zeros(49,Con2); 

for i=1:Con2
    if i==1
        
        par.Pcr=0;
        par.Ptr=0;
        
    elseif i==2
        
        par.Pcr=0.1*10^(-9);
        par.Ptr=0.1*10^(-9);
       
    else
    
        par.Pcr=0.25*10^(-9);
        par.Ptr=0.25*10^(-9);
    

    end
    tspan=0:300:14400;   
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    y0=[0 0 0 par.dCas9tot 0 0 0 0 0 0];  
    [t,y]=ode23s(@CRISPRModel,tspan,y0,options,par);
    y=y*10^6;
    Simu(:,i)=y(:,10);
    
end

OError=sum((Simu(:,1)-Ave09(:,1)).^2)+sum((Simu(:,2)-Ave09(:,2)).^2)...
    +sum((Simu(:,3)-Ave09(:,3)).^2);



iter=1;
Maxiter=210000; 
en=0;
TE=zeros(Maxiter,1);
Rpara=zeros(Maxiter,12);


while iter<Maxiter  
    clear y Newpar
    Newpar=par;
    step=randi([-1 1],1,12);
    Newpar.alpha_cr=par.alpha_cr+step(1)*stepRalpha;
    Newpar.alpha_m=par.alpha_m+step(2)*stepRaM;
    Newpar.deg_tr=par.deg_tr+step(3)*stepDeg;
    Newpar.deg_cr=par.deg_cr+step(4)*stepDeg;
    Newpar.deg_g=par.deg_g+step(5)*stepDeg;
    Newpar.deg_m=par.deg_m+step(6)*stepDeg;
    Newpar.gamma1=par.gamma1+step(7)*stepGamma;
    Newpar.gamma2=par.gamma2+step(8)*stepGamma;
    Newpar.omega=par.omega+step(9)*stepOmega;  
    Newpar.KI=par.KI+step(10)*stepKi; 
    Newpar.KE=par.KE+step(11)*stepKe; 
    Newpar.alpha_gm=par.alpha_gm+step(12)*stepGm; 


    if (Newpar.alpha_cr>=Lalpha) && (Newpar.alpha_cr<=Ualpha)...
        && (Newpar.alpha_m>=LaM) && (Newpar.alpha_m<=UaM)...
        && (Newpar.gamma1>=Lgamma) && (Newpar.gamma1<=Ugamma)...
        && (Newpar.gamma2>=Lgamma) && (Newpar.gamma2<=Ugamma)...
        && (Newpar.deg_tr>=Ldeg) && (Newpar.deg_tr<=Udeg)...
        && (Newpar.deg_cr>=Ldeg) && (Newpar.deg_cr<=Udeg)...
        && (Newpar.deg_g>=Ldeg) && (Newpar.deg_g<=Udeg)...
        && (Newpar.deg_m>=Ldeg) && (Newpar.deg_m<=Udeg)...
        && (Newpar.omega>=Lomega) && (Newpar.omega<=Uomega)... 
        && (Newpar.KI>=LKi) && (Newpar.KI<=UKi)...
        && (Newpar.KE>=LKe) && (Newpar.KE<=UKe)...
        && (Newpar.alpha_gm>=LGm) && (Newpar.alpha_gm<=UGm)
    
    
    en=en+1;
    Simu=zeros(49,Con2);    
        for i=1:Con2
            if i==1
                
                Newpar.Pcr=0;
                Newpar.Ptr=0;

            elseif i==2
                
                Newpar.Pcr=0.1*10^(-9);
                Newpar.Ptr=0.1*10^(-9);
                
            else

           
                Newpar.Pcr=0.25*10^(-9);
                Newpar.Ptr=0.25*10^(-9);

 
            end
            tspan=0:300:14400;  
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            y0=[0 0 0 par.dCas9tot 0 0 0 0 0 0];  
            [t,y]=ode23s(@CRISPRModel,tspan,y0,options,Newpar);
            y=y*10^6;
            Simu(:,i)=y(:,10);
            
        end
        NError=sum((Simu(:,1)-Ave09(:,1)).^2)+sum((Simu(:,2)-Ave09(:,2)).^2)...
                +sum((Simu(:,3)-Ave09(:,3)).^2);
           
        Pro=exp(-(NError-OError)/(2*0.25^2));  
        
        if NError<OError
            OError=NError;
            par=Newpar;
            iter=iter+1;
            TE(iter)=OError;
        elseif Pro>rand(1)
            OError=NError;
            par=Newpar;
            iter=iter+1;
            TE(iter)=OError;
        else
            iter=iter+1;
            TE(iter)=OError;
        end
        
    else
    iter=iter+1;
    TE(iter)=OError;
        
    end
    
    Rpara(iter,1)=par.alpha_gm;
    Rpara(iter,2)=par.alpha_cr;
    Rpara(iter,3)=par.alpha_m;
    Rpara(iter,4)=par.omega;
    Rpara(iter,5)=par.gamma1;
    Rpara(iter,6)=par.gamma2;
    Rpara(iter,7)=par.deg_tr;
    Rpara(iter,8)=par.deg_cr;
    Rpara(iter,9)=par.deg_g;
    Rpara(iter,10)=par.deg_m;
    Rpara(iter,11)=par.KI;
    Rpara(iter,12)=par.KE;
    
    
  clear y  
end

[v,loca]=min(TE(2:end));
par.alpha_gm=Rpara(loca+1,1);
par.alpha_cr=Rpara(loca+1,2);
par.alpha_m=Rpara(loca+1,3);
par.omega=Rpara(loca+1,4);
par.gamma1=Rpara(loca+1,5);
par.gamma2=Rpara(loca+1,6);
par.deg_tr=Rpara(loca+1,7);
par.deg_cr=Rpara(loca+1,8);
par.deg_g=Rpara(loca+1,9);
par.deg_m=Rpara(loca+1,10);
par.KI=Rpara(loca+1,11);
par.KE=Rpara(loca+1,12);


Simu=zeros(49,Con2); 

for i=1:Con2
    if i==1
        
        par.Pcr=0;
        par.Ptr=0;
        
    elseif i==2
        
        par.Pcr=0.1*10^(-9);
        par.Ptr=0.1*10^(-9);
    else   

        par.Pcr=0.25*10^(-9);
        par.Ptr=0.25*10^(-9);
    
    end
    tspan=0:300:14400;   
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    y0=[0 0 0 par.dCas9tot 0 0 0 0 0 0];  
    [t,y]=ode23s(@CRISPRModel,tspan,y0,options,par);
    y=y*10^6;
    Simu(:,i)=y(:,10);
    
end

save(savefile,'Rpara','par','TE','en')

subplot(3,4,nj)
plot(tspan./60,Simu(:,1),'-r',tspan./60,Ave09(:,1),'-.r',...
    tspan./60,Simu(:,2),'-b',tspan./60,Ave09(:,2),'-.b',...
    tspan./60,Simu(:,3),'-y',tspan./60,Ave09(:,3),'-.y','LineWidth',2)
xlabel('Time (min)')
ylabel('Concentration (\muM)')
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')


end
