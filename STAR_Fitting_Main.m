%%%% simulates the titration process


close all
clear all

Con2=3;

M11=csvread('Data3011.csv',14,0);
Ave11=zeros(49,Con2);
for i=1:Con2
    for j=1:9
        Ave11(:,i)=Ave11(:,i)+M11(:,(i-1)*9+j+1);
    end
end
Ave11=Ave11./9;


%%% parameter intervals

Lalpha=10^(-2);
Ualpha=10;
stepRalpha=(Ualpha-Lalpha)/1000;  %%% for transcription

Ldeg=10^(-5);
Udeg=10^(-1);
stepDeg=(Udeg-Ldeg)/1000; %%% for degradation

Lbd=10^3;
Ubd=10^7;
stepBeta=(Ubd-Lbd)/1000; %%% for STAR binding

LaM=10^(-2); %%% original 10^-1
UaM=100;
stepRaM=(UaM-LaM)/1000; %%% GFP mRNA transcription rate

LKi=10^(-4);
UKi=10^(-2);
stepKi=(UKi-LKi)/1000; %%% initiation rate

LKe=10^(-4);
UKe=10^(-2);
stepKe=(UKe-LKe)/1000; %%% elongation rate

LGm=10^(-3);
UGm=10^(-1);
stepGm=(UGm-LGm)/1000; %%% maturation rate

%%%% fitting 8 parameters par.alpha_s; par.deg_s; par.beta_s
%%%% par.alpha_m; par.deg_m; par.KI; par.KE; par.alpha_gm


figure
Percent=[0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95];

for nj=1:10
    savefile=strcat('M11Fit105000Final',int2str(nj),'.mat');
    
    Sample=randi([1 10],1,8);
       
    par.alpha_s=Lalpha+(Ualpha-Lalpha)*Percent(Sample(1));
    par.alpha_m=LaM+(UaM-LaM)*Percent(Sample(2));
    par.beta_s=Lbd+(Ubd-Lbd)*Percent(Sample(3));
    par.deg_s=Ldeg+(Udeg-Ldeg)*Percent(Sample(4));
    par.deg_m=Ldeg+(Udeg-Ldeg)*Percent(Sample(5));
    par.KI=LKi+(UKi-LKi)*Percent(Sample(6)); %%% initiation . 10^(-4) to 10^(-2)
    par.KE=LKe+(UKe-LKe)*Percent(Sample(7)); %%% elongation . 10^(-4) to 10^(-2)
    par.alpha_gm=LGm+(UGm-LGm)*Percent(Sample(8)); %%% maturation . 10^(-3) to 10^(-1)

   
    
par.Pytot=0.5*10^(-9);  %%% 0.5nM STAR targeting GFP plasmid


%%%%

Simu=zeros(49,Con2); %%% every 5 mins measurement

for i=2:Con2
    
    if i==1  %%% control experiments
        par.Ps=0;
        
    elseif i==2
        par.Ps=4*10^(-9); %%% STAR
        
    else
        par.Ps=8*10^(-9); %%% STAR
        
    end
    
    tspan=0:300:14400; %%% 5 mins increment
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    y0=[0 0 0 0 0];  
    [t,y]=ode23s(@STARModel,tspan,y0,options,par); %%% call the model that has activation 
    y=y*10^6;
    
Simu(:,i)=y(:,5);
end
OError=sum((Simu(:,2)-Ave11(:,2)).^2)...
    +sum((Simu(:,3)-Ave11(:,3)).^2);




iter=1;
Maxiter=105000;  %%% max iteration
en=0;

TE=zeros(Maxiter,1);
Rpara=zeros(Maxiter,8);

while iter<Maxiter
    
    clear y Newpar
    Newpar=par;
    
    step=randi([-1 1],1,8);
    Newpar.alpha_s=par.alpha_s+step(1)*stepRalpha;
    Newpar.deg_s=par.deg_s+step(2)*stepDeg;
    Newpar.beta_s=par.beta_s+step(3)*stepBeta; 
    Newpar.alpha_m=par.alpha_m+step(4)*stepRaM;
    Newpar.deg_m=par.deg_m+step(5)*stepDeg;
    Newpar.KI=par.KI+step(6)*stepKi; 
    Newpar.KE=par.KE+step(7)*stepKe; 
    Newpar.alpha_gm=par.alpha_gm+step(8)*stepGm; %%% maturation . 10^(-5) to 10^(-3)
        
    if (Newpar.alpha_s>=Lalpha) && (Newpar.alpha_s<=Ualpha)...
        && (Newpar.alpha_m>=LaM) && (Newpar.alpha_m<=UaM)...
        && (Newpar.deg_s>=Ldeg) && (Newpar.deg_s<=Udeg)...
        && (Newpar.beta_s>=Lbd) && (Newpar.beta_s<=Ubd)...
        && (Newpar.deg_m>=Ldeg) && (Newpar.deg_m<=Udeg)...
        && (Newpar.KI>=LKi) && (Newpar.KI<=UKi)...
        && (Newpar.KE>=LKe) && (Newpar.KE<=UKe)...
        && (Newpar.alpha_gm>=LGm) && (Newpar.alpha_gm<=UGm)
        
        en=en+1; 
        Simu=zeros(49,Con2); %%% every 5 mins measurement
        for i=2:Con2
            if i==1  %%% control experiments
                Newpar.Ps=0;
                
            elseif i==2
                Newpar.Ps=4*10^(-9); %%% STAR
                
            else
                Newpar.Ps=8*10^(-9); %%% STAR        
            end

            tspan=0:300:14400; %%% 5 mins increment
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            y0=[0 0 0 0 0];  
            [t,y]=ode23s(@STARModel,tspan,y0,options,Newpar);
            y=y*10^6;
            Simu(:,i)=y(:,5);
        end

        NError=sum((Simu(:,2)-Ave11(:,2)).^2)...
               +sum((Simu(:,3)-Ave11(:,3)).^2);
        Pro=exp(-(NError-OError)/(2*0.25^2));  %%% Metropolis criteria with T=0.125 
                                               %%% with standard deviation as 0.25
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
        
    Rpara(iter,1)=par.alpha_s;
    Rpara(iter,2)=par.alpha_m;
    Rpara(iter,3)=par.deg_s;
    Rpara(iter,4)=par.beta_s;
    Rpara(iter,5)=par.deg_m;
    Rpara(iter,6)=par.KI;
    Rpara(iter,7)=par.KE;
    Rpara(iter,8)=par.alpha_gm;
    
  clear y  
end

[v,loca]=min(TE(2:end));

par.alpha_s=Rpara(loca+1,1);
par.alpha_m=Rpara(loca+1,2);
par.deg_s=Rpara(loca+1,3);
par.beta_s=Rpara(loca+1,4);
par.deg_m=Rpara(loca+1,5);
par.KI=Rpara(loca+1,6);
par.KE=Rpara(loca+1,7);
par.alpha_gm=Rpara(loca+1,8);

Simu=zeros(49,Con2); %%% every 5 mins measurement
for i=2:Con2
    if i==1  %%% control experiments
        par.Ps=0;
        
    elseif i==2
        par.Ps=4*10^(-9); %%% STAR
       
    else
        par.Ps=8*10^(-9); %%% STAR
        
    end

    tspan=0:300:14400;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    y0=[0 0 0 0 0];  
    [t,y]=ode23s(@STARModel,tspan,y0,options,par);
    y=y*10^6;
    
Simu(:,i)=y(:,5);
end


save(savefile,'Rpara','par','en','TE')

subplot(3,4,nj)
plot(tspan./60,Simu(:,1),'-r',tspan./60,Ave11(:,1),'-.r',...
    tspan./60,Simu(:,2),'-b',tspan./60,Ave11(:,2),'-.b',...
    tspan./60,Simu(:,3),'-m',tspan./60,Ave11(:,3),'-.m','LineWidth',2)
xlabel('Time (min)')
ylabel('Concentration (\muM)')
xlim([0 250])
ylim([0 0.1])
% % % legend('SNT','ENT','SRpr','ERpr')
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

% % % subplot(4,2,nj*2)
% % % plot(1:iter,TE(1:iter),'-k','LineWidth',2)
% % % xlabel('Iteration')
% % % ylabel('Error (\muM)')
% % % set(gca,'FontSize',24)
% % % set(gca,'FontName','Times New Roman')

    
end
    


