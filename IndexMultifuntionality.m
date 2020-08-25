%% Code used to plot the figures found in the article "An infrastructure perspective for\\
%% enhancing multi-functionality of forests:A conceptual modeling approach" submitted to the Earths'
%% The code simulates the proposed model in the article according to the defined multi-functionality index (Figure 6 in the manuescript)


% clc;
% clear all;
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=0.75;
s=0.0125;
u=0.0167;
h=0.025;
g1=0.16;d1=0.45;h1=25;v1=2.29;
g2=0.013;d2=0.125;h2=7;v2=0.066;
z=0.0008;
d=0.0067;
alpha=0.06;
r=0.015;
pe=0.9;
pa=0.9;
pd=0.8;
p=50;
cm=20;
gamma=0.1;

phi=0.2;
SIm=0.9;
SI0=0.1;
SIm1=0.4;
SI01=0.01;
alpha1=0.5;
u1=0.2*20;
alpha2=0.2;
u2=0.2;
beta=0.1;
k=80;
delta=0.05;
alphaT=0.00009;
alphaC=0.01;
CI0=0.5;
aT=1;
CT=1;
TCF=0.4;
TCT=TCF*0.5;
pde=40;
cde=20;

pyT=1;
j=25;
counter=0;
t1=[];
t2=[];
L=[];
Y=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oo=2;

U=1;
Q=1;
mu=0.0005;
kk=10;
II=0.5;

HH=@(x,y,SI0,SIm) (0).*(x<SI0) + (y.*((x-SI0)/(SIm-SI0))).*(x>SI0 & x<SIm) + (y).*(x>SIm);
HHT=@(a,b,c,d,e,f,g) (0).*(b<e) + (f*exp((-1/2)*((a-c)/g)^2).*((b-e)/(d-e))).*(b>e & b<d) + (f*exp((-1/2)*((a-c)/g)^2)).*(b>d);
% HHT(y4(t-1),y5(t-1),CI0,SIm,SI0,aT,CT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ooo=[];
iii=[];
NN=[];
PP=[];
AM=[];
XX=1;
AC=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ee=15;
dt=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation for a tax-investment tradeoff (Mean && Diff)
for amain=1:60
    amainn=0.01+(amain-1)/150;
    amain
    for acons=1:60
        aconss=0.02+(acons-1)/60;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y1(1)=Q*200;y2(1)=Q*300;y3(1)=25*Q;y4(1)=0.3;y5(1)=0.4;y6(1)=1000;y7(1)=0.1;BA(1)=0.07;E(1)=0.01;P(1)=0.001;II=0.5;hh=kk*Q;o=oo*Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ooo=[ooo;amainn];
            iii=[iii;aconss];
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
TCT=amainn*0.5;
TCF=amainn;
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
wood=[];
BB=[];
minstrata1=[];

                 for t=2:500                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forest dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Strata 1
                        if y1(t-1)>0
                            hh=kk*Q;
                            y1(t) = y1(t-1) + dt*(h*y2(t-1)*(1-u*g1*y1(t-1)) - y1(t-1)*d - hh*HH(y5(t-1),y4(t-1),SI0,SIm)/(v1));
                        else
                            y1(t)=0 ;
                            hh=0;
                            II=0;
                            XX=0;
                        end 
                %%% Strata 2
                        if y2(t-1)>0
                            y2(t) = y2(t-1) + dt*(b*g1*y1(t-1)*(1-s*(g1*y1(t-1) + g2*y2(t-1))) - h*y2(t-1)*(1-u*g1*y1(t-1)) - y2(t-1)*(z*g1*y1(t-1)+d));
                        else
                            y2(t)=0;
                        end
                %%% deadwood dynamics
                        if y3(t-1)>0
                            o=oo*Q;
                            y3(t) = y3(t-1) + dt*(v2*y2(t-1)*(z*g1*y1(t-1)+d) + (hh*HH(y5(t-1),y4(t-1),SI0,SIm))*(1-pe) - (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm)) ...
                             + v1*y1(t-1)*d - alpha*y3(t-1));
                        else
                            y3(t) = 0;
                            o=0;
                        end
                %%% Annual budget
                        BA(t)= dt*(U*(TCF*(hh*(HH(y5(t-1),y4(t-1),SI0,SIm))*(p-cm)+...
                                   (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))*(pde-cde))+gamma+TCT*y6(t-1)));                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Infrastructures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Infrastructure construction dynamics
                        if y4(t-1)>0
                            y4(t) = y4(t-1) +   dt*(II*(aconss)*BA(t)*u1)*mu*y4(t-1)*(1-(y4(t-1)/1)); 
                        else
                            y4(t)=0;
                        end 

                %%% Infrastructures state dynamics
                            ez=((1-aconss)*II*BA(t)*u2*(1/(y4(t-1)+k)));
                            y5(t) = y5(t-1) + dt*((ez)*y5(t-1)*((1-y5(t-1))/1)+(y4(t)-y4(t-1))*((1-y5(t-1))/y4(t-1))- delta*y5(t-1)) ; %%logistic equation          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tourism dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sigma1 = 1;
                sigma2 = 1;
                scale1 = 100;
                scale2 = 50;
                sigma1 = scale1*sigma1;
                sigma2 = scale2*sigma2;
                Theta = -55;

                aaa = ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
                bbb = -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
                ccc = ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));

                muu = [200 400];
                AAA = 1;
                E = AAA*exp(-(aaa*(y1(t-1) - muu(1)).^2 + 2*bbb*(y1(t-1) - muu(1)).*(y2(t-1) - muu(2)) + ccc*(y2(t-1) - muu(2)).^2));
                
                        if y6(t-1)>0
                            if y1(t-1)>20
                                alphaC=0.1;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            else
                                alphaC=0.15;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            end
                        else
                            y6(t)=0;
                        end

                    EE=[EE,E];
                    PP=[PP,HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT)];
                    CCF=[CCF,HH(y5(t-1),y4(t-1),SI0,SIm)];
                    CC=[CC,-alphaT*y6(t-1)];
                    TT=[TT,t];
                    MM=[MM,II*(alpha1*BA(t)*u1*(1/(y4(t-1)+k)))];
                    MMM=[MMM,-delta*y5(t-1)];
                    CCT=[CCT,(aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2))];
                    wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
                    BB=[BB,y1(t)-y1(t-1)];
                 end
                 wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
               
%% Diff
            Strata1D1(amain,acons) = y1(500);            
            Strata2D1(amain,acons) = y2(500);
            DeadD1(amain,acons) = y3(500);
            COND1(amain,acons) = y4(500);
            STATED1(amain,acons) = y5(500);
            TOURISMD1(amain,acons) = y6(500);
            BUDGETD1(amain,acons) = BA(500);
            ENHANCD1(amain,acons) = CCF(499);
            ENVID1(amain,acons) = EE(499);
            WOODD1(amain,acons) = wood(499);
            AINFRAD1(amain,acons) = PP(499);

    end


end


            Strata1D1 = Strata1D1/max(Strata1D1(:));            
            Strata2D1 = Strata2D1/max(Strata2D1(:));
            DeadD1 = DeadD1/max(DeadD1(:));
            COND1 = COND1/max(COND1(:));
            STATED1 = STATED1/max(STATED1(:));
            TOURISMD1 = TOURISMD1/max(TOURISMD1(:));
            BUDGETD1 = BUDGETD1/max(BUDGETD1(:));
            ENHANCD1 = ENHANCD1/max(ENHANCD1(:));
            ENVID1 = ENVID1/max(ENVID1(:));
            WOODD1 = WOODD1/max(WOODD1(:));
            AINFRAD1 =  AINFRAD1/max(AINFRAD1(:));  
            

            
                        %% multifunctionality indx with government
for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex1(amainl,aconsl)=0.333*WOODD1(amainl,aconsl)+0.333*TOURISMD1(amainl,aconsl)+0.333*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex2(amainl,aconsl)=0.5*WOODD1(amainl,aconsl)+0.25*TOURISMD1(amainl,aconsl)+0.25*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex3(amainl,aconsl)=0.25*WOODD1(amainl,aconsl)+0.5*TOURISMD1(amainl,aconsl)+0.25*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex4(amainl,aconsl)=0.25*WOODD1(amainl,aconsl)+0.25*TOURISMD1(amainl,aconsl)+0.5*DeadD1(amainl,aconsl);   
    end
end



    %% function management foster   %%%%%%%%%% figure 4
figure(1)
colormap(jet)
grid on
hold on
subplot(2,4,1)
imagesc(Mindex1)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,2)
imagesc(Mindex2)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,3)
imagesc(Mindex3)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,4)
imagesc(Mindex4)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=0.75;
s=0.0125;
u=0.0167;
h=0.025;
g1=0.16;d1=0.45;h1=25;v1=2.29;
g2=0.013;d2=0.125;h2=7;v2=0.066;
z=0.0008;
d=0.0067;
alpha=0.06;
r=0.015;
pe=0.9;
pa=0.9;
pd=0.8;
p=50;
cm=20;
gamma=0.1;

phi=0.2;
SIm=0.9;
SI0=0.1;
SIm1=0.4;
SI01=0.01;
alpha1=0.5;
u1=0.2*20;
alpha2=0.2;
u2=0.2;
beta=0.1;
k=80;
delta=0.05;
alphaT=0.00009;
alphaC=0.01;
CI0=0.5;
aT=1;
CT=1;
TCF=0.4;
TCT=TCF*0.5;
pde=40;
cde=20;

pyT=1;
j=25;
counter=0;
t1=[];
t2=[];
L=[];
Y=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oo=3;

U=1;
Q=1;
mu=0.0005;
kk=30;
II=0.5;

HH=@(x,y,SI0,SIm) (0).*(x<SI0) + (y.*((x-SI0)/(SIm-SI0))).*(x>SI0 & x<SIm) + (y).*(x>SIm);
HHT=@(a,b,c,d,e,f,g) (0).*(b<e) + (f*exp((-1/2)*((a-c)/g)^2).*((b-e)/(d-e))).*(b>e & b<d) + (f*exp((-1/2)*((a-c)/g)^2)).*(b>d);
% HHT(y4(t-1),y5(t-1),CI0,SIm,SI0,aT,CT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ooo=[];
iii=[];
NN=[];
PP=[];
AM=[];
XX=1;
AC=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ee=15;
dt=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation for a tax-investment tradeoff (Mean && Diff)
for amain=1:60
    amainn=0.01+(amain-1)/150;
    amain
    for acons=1:60
        aconss=0.02+(acons-1)/60;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y1(1)=Q*200;y2(1)=Q*300;y3(1)=25*Q;y4(1)=0.3;y5(1)=0.4;y6(1)=1000;y7(1)=0.1;BA(1)=0.07;E(1)=0.01;P(1)=0.001;II=0.5;hh=kk*Q;o=oo*Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ooo=[ooo;amainn];
            iii=[iii;aconss];
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
TCT=amainn*0.5;
TCF=amainn;
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
wood=[];
BB=[];
minstrata1=[];

                 for t=2:500                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forest dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Strata 1
                        if y1(t-1)>0
                            hh=kk*Q;
                            y1(t) = y1(t-1) + dt*(h*y2(t-1)*(1-u*g1*y1(t-1)) - y1(t-1)*d - hh*HH(y5(t-1),y4(t-1),SI0,SIm)/(v1));
                        else
                            y1(t)=0 ;
                            hh=0;
                            II=0;
                            XX=0;
                        end 
                %%% Strata 2
                        if y2(t-1)>0
                            y2(t) = y2(t-1) + dt*(b*g1*y1(t-1)*(1-s*(g1*y1(t-1) + g2*y2(t-1))) - h*y2(t-1)*(1-u*g1*y1(t-1)) - y2(t-1)*(z*g1*y1(t-1)+d));
                        else
                            y2(t)=0;
                        end
                %%% deadwood dynamics
                        if y3(t-1)>0
                            o=oo*Q;
                            y3(t) = y3(t-1) + dt*(v2*y2(t-1)*(z*g1*y1(t-1)+d) + (hh*HH(y5(t-1),y4(t-1),SI0,SIm))*(1-pe) - (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm)) ...
                             + v1*y1(t-1)*d - alpha*y3(t-1));
                        else
                            y3(t) = 0;
                            o=0;
                        end
                %%% Annual budget
                        BA(t)= dt*(U*(TCF*(hh*(HH(y5(t-1),y4(t-1),SI0,SIm))*(p-cm)+...
                                   (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))*(pde-cde))+gamma+TCT*y6(t-1)));                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Infrastructures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Infrastructure construction dynamics
                        if y4(t-1)>0
                            y4(t) = y4(t-1) +   dt*(II*(aconss)*BA(t)*u1)*mu*y4(t-1)*(1-(y4(t-1)/1)); 
                        else
                            y4(t)=0;
                        end 

                %%% Infrastructures state dynamics
                            ez=((1-aconss)*II*BA(t)*u2*(1/(y4(t-1)+k)));
                            y5(t) = y5(t-1) + dt*((ez)*y5(t-1)*((1-y5(t-1))/1)+(y4(t)-y4(t-1))*((1-y5(t-1))/y4(t-1))- delta*y5(t-1)) ; %%logistic equation          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tourism dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sigma1 = 1;
                sigma2 = 1;
                scale1 = 100;
                scale2 = 50;
                sigma1 = scale1*sigma1;
                sigma2 = scale2*sigma2;
                Theta = -55;

                aaa = ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
                bbb = -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
                ccc = ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));

                muu = [200 400];
                AAA = 1;
                E = AAA*exp(-(aaa*(y1(t-1) - muu(1)).^2 + 2*bbb*(y1(t-1) - muu(1)).*(y2(t-1) - muu(2)) + ccc*(y2(t-1) - muu(2)).^2));
                
                        if y6(t-1)>0
                            if y1(t-1)>20
                                alphaC=0.1;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            else
                                alphaC=0.15;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            end
                        else
                            y6(t)=0;
                        end

                    EE=[EE,E];
                    PP=[PP,HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT)];
                    CCF=[CCF,HH(y5(t-1),y4(t-1),SI0,SIm)];
                    CC=[CC,-alphaT*y6(t-1)];
                    TT=[TT,t];
                    MM=[MM,II*(alpha1*BA(t)*u1*(1/(y4(t-1)+k)))];
                    MMM=[MMM,-delta*y5(t-1)];
                    CCT=[CCT,(aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2))];
                    wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
                    BB=[BB,y1(t)-y1(t-1)];
                 end
                 wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
               
%% Diff
            Strata1D1(amain,acons) = y1(500);            
            Strata2D1(amain,acons) = y2(500);
            DeadD1(amain,acons) = y3(500);
            COND1(amain,acons) = y4(500);
            STATED1(amain,acons) = y5(500);
            TOURISMD1(amain,acons) = y6(500);
            BUDGETD1(amain,acons) = BA(500);
            ENHANCD1(amain,acons) = CCF(499);
            ENVID1(amain,acons) = EE(499);
            WOODD1(amain,acons) = wood(499);
            AINFRAD1(amain,acons) = PP(499);

    end


end


            Strata1D1 = Strata1D1/max(Strata1D1(:));            
            Strata2D1 = Strata2D1/max(Strata2D1(:));
            DeadD1 = DeadD1/max(DeadD1(:));
            COND1 = COND1/max(COND1(:));
            STATED1 = STATED1/max(STATED1(:));
            TOURISMD1 = TOURISMD1/max(TOURISMD1(:));
            BUDGETD1 = BUDGETD1/max(BUDGETD1(:));
            ENHANCD1 = ENHANCD1/max(ENHANCD1(:));
            ENVID1 = ENVID1/max(ENVID1(:));
            WOODD1 = WOODD1/max(WOODD1(:));
            AINFRAD1 =  AINFRAD1/max(AINFRAD1(:));  
            

            
                        %% multifunctionality indx with government
for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex1(amainl,aconsl)=0.333*WOODD1(amainl,aconsl)+0.333*TOURISMD1(amainl,aconsl)+0.333*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex2(amainl,aconsl)=0.5*WOODD1(amainl,aconsl)+0.25*TOURISMD1(amainl,aconsl)+0.25*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex3(amainl,aconsl)=0.25*WOODD1(amainl,aconsl)+0.5*TOURISMD1(amainl,aconsl)+0.25*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex4(amainl,aconsl)=0.25*WOODD1(amainl,aconsl)+0.25*TOURISMD1(amainl,aconsl)+0.5*DeadD1(amainl,aconsl);   
    end
end



    %% function management foster   %%%%%%%%%% figure 4
figure(2)
colormap(jet)
grid on
hold on
subplot(2,4,1)
imagesc(Mindex1)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,2)
imagesc(Mindex2)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,3)
imagesc(Mindex3)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,4)
imagesc(Mindex4)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=0.75;
s=0.0125;
u=0.0167;
h=0.025;
g1=0.16;d1=0.45;h1=25;v1=2.29;
g2=0.013;d2=0.125;h2=7;v2=0.066;
z=0.0008;
d=0.0067;
alpha=0.06;
r=0.015;
pe=0.9;
pa=0.9;
pd=0.8;
p=50;
cm=20;
gamma=0.1;

phi=0.2;
SIm=0.9;
SI0=0.1;
SIm1=0.4;
SI01=0.01;
alpha1=0.5;
u1=0.2*20;
alpha2=0.2;
u2=0.2;
beta=0.1;
k=80;
delta=0.05;
alphaT=0.00005;
alphaC=0.01;
CI0=0.5;
aT=1;
CT=1;
TCF=0.4;
TCT=TCF*0.5;
pde=40;
cde=20;

pyT=1;
j=25;
counter=0;
t1=[];
t2=[];
L=[];
Y=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oo=2;

U=1;
Q=1;
mu=0.0005;
kk=10;
II=0.5;

HH=@(x,y,SI0,SIm) (0).*(x<SI0) + (y.*((x-SI0)/(SIm-SI0))).*(x>SI0 & x<SIm) + (y).*(x>SIm);
HHT=@(a,b,c,d,e,f,g) (0).*(b<e) + (f*exp((-1/2)*((a-c)/g)^2).*((b-e)/(d-e))).*(b>e & b<d) + (f*exp((-1/2)*((a-c)/g)^2)).*(b>d);
% HHT(y4(t-1),y5(t-1),CI0,SIm,SI0,aT,CT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ooo=[];
iii=[];
NN=[];
PP=[];
AM=[];
XX=1;
AC=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ee=15;
dt=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation for a tax-investment tradeoff (Mean && Diff)
for amain=1:60
    amainn=0.01+(amain-1)/150;
    amain
    for acons=1:60
        aconss=0.02+(acons-1)/60;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y1(1)=Q*200;y2(1)=Q*300;y3(1)=25*Q;y4(1)=0.3;y5(1)=0.4;y6(1)=1000;y7(1)=0.1;BA(1)=0.07;E(1)=0.01;P(1)=0.001;II=0.5;hh=kk*Q;o=oo*Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ooo=[ooo;amainn];
            iii=[iii;aconss];
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
TCT=amainn*0.5;
TCF=amainn;
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
wood=[];
BB=[];
minstrata1=[];

                 for t=2:500                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forest dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Strata 1
                        if y1(t-1)>0
                            hh=kk*Q;
                            y1(t) = y1(t-1) + dt*(h*y2(t-1)*(1-u*g1*y1(t-1)) - y1(t-1)*d - hh*HH(y5(t-1),y4(t-1),SI0,SIm)/(v1));
                        else
                            y1(t)=0 ;
                            hh=0;
                            II=0;
                            XX=0;
                        end 
                %%% Strata 2
                        if y2(t-1)>0
                            y2(t) = y2(t-1) + dt*(b*g1*y1(t-1)*(1-s*(g1*y1(t-1) + g2*y2(t-1))) - h*y2(t-1)*(1-u*g1*y1(t-1)) - y2(t-1)*(z*g1*y1(t-1)+d));
                        else
                            y2(t)=0;
                        end
                %%% deadwood dynamics
                        if y3(t-1)>0
                            o=oo*Q;
                            y3(t) = y3(t-1) + dt*(v2*y2(t-1)*(z*g1*y1(t-1)+d) + (hh*HH(y5(t-1),y4(t-1),SI0,SIm))*(1-pe) - (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm)) ...
                             + v1*y1(t-1)*d - alpha*y3(t-1));
                        else
                            y3(t) = 0;
                            o=0;
                        end
                %%% Annual budget
                        BA(t)= dt*(U*(TCF*(hh*(HH(y5(t-1),y4(t-1),SI0,SIm))*(p-cm)+...
                                   (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))*(pde-cde))+gamma+TCT*y6(t-1)));                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Infrastructures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Infrastructure construction dynamics
                        if y4(t-1)>0
                            y4(t) = y4(t-1) +   dt*(II*(aconss)*BA(t)*u1)*mu*y4(t-1)*(1-(y4(t-1)/1)); 
                        else
                            y4(t)=0;
                        end 

                %%% Infrastructures state dynamics
                            ez=((1-aconss)*II*BA(t)*u2*(1/(y4(t-1)+k)));
                            y5(t) = y5(t-1) + dt*((ez)*y5(t-1)*((1-y5(t-1))/1)+(y4(t)-y4(t-1))*((1-y5(t-1))/y4(t-1))- delta*y5(t-1)) ; %%logistic equation          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tourism dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sigma1 = 1;
                sigma2 = 1;
                scale1 = 100;
                scale2 = 50;
                sigma1 = scale1*sigma1;
                sigma2 = scale2*sigma2;
                Theta = -55;

                aaa = ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
                bbb = -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
                ccc = ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));

                muu = [200 400];
                AAA = 1;
                E = AAA*exp(-(aaa*(y1(t-1) - muu(1)).^2 + 2*bbb*(y1(t-1) - muu(1)).*(y2(t-1) - muu(2)) + ccc*(y2(t-1) - muu(2)).^2));
                
                        if y6(t-1)>0
                            if y1(t-1)>20
                                alphaC=0.1;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            else
                                alphaC=0.15;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            end
                        else
                            y6(t)=0;
                        end

                    EE=[EE,E];
                    PP=[PP,HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT)];
                    CCF=[CCF,HH(y5(t-1),y4(t-1),SI0,SIm)];
                    CC=[CC,-alphaT*y6(t-1)];
                    TT=[TT,t];
                    MM=[MM,II*(alpha1*BA(t)*u1*(1/(y4(t-1)+k)))];
                    MMM=[MMM,-delta*y5(t-1)];
                    CCT=[CCT,(aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2))];
                    wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
                    BB=[BB,y1(t)-y1(t-1)];
                 end
                 wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
               
%% Diff
            Strata1D1(amain,acons) = y1(500);            
            Strata2D1(amain,acons) = y2(500);
            DeadD1(amain,acons) = y3(500);
            COND1(amain,acons) = y4(500);
            STATED1(amain,acons) = y5(500);
            TOURISMD1(amain,acons) = y6(500);
            BUDGETD1(amain,acons) = BA(500);
            ENHANCD1(amain,acons) = CCF(499);
            ENVID1(amain,acons) = EE(499);
            WOODD1(amain,acons) = wood(499);
            AINFRAD1(amain,acons) = PP(499);

    end


end


            Strata1D1 = Strata1D1/max(Strata1D1(:));            
            Strata2D1 = Strata2D1/max(Strata2D1(:));
            DeadD1 = DeadD1/max(DeadD1(:));
            COND1 = COND1/max(COND1(:));
            STATED1 = STATED1/max(STATED1(:));
            TOURISMD1 = TOURISMD1/max(TOURISMD1(:));
            BUDGETD1 = BUDGETD1/max(BUDGETD1(:));
            ENHANCD1 = ENHANCD1/max(ENHANCD1(:));
            ENVID1 = ENVID1/max(ENVID1(:));
            WOODD1 = WOODD1/max(WOODD1(:));
            AINFRAD1 =  AINFRAD1/max(AINFRAD1(:));  
            

            
             %% multifunctionality indx with government
             
             
for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex1(amainl,aconsl)=0.333*WOODD1(amainl,aconsl)+0.333*TOURISMD1(amainl,aconsl)+0.333*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex2(amainl,aconsl)=0.5*WOODD1(amainl,aconsl)+0.25*TOURISMD1(amainl,aconsl)+0.25*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex3(amainl,aconsl)=0.25*WOODD1(amainl,aconsl)+0.5*TOURISMD1(amainl,aconsl)+0.25*DeadD1(amainl,aconsl);   
    end
end

for amainl=1:60
    amainl
    AM=[AM;amainn];
    for aconsl=1:60
        AC=[AC;aconss];
       Mindex4(amainl,aconsl)=0.25*WOODD1(amainl,aconsl)+0.25*TOURISMD1(amainl,aconsl)+0.5*DeadD1(amainl,aconsl);   
    end
end



    %% function management foster   %%%%%%%%%% figure 4
figure(3)
colormap(jet)
grid on
hold on
subplot(2,4,1)
imagesc(Mindex1)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,2)
imagesc(Mindex2)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,3)
imagesc(Mindex3)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on
hold on
subplot(2,4,4)
imagesc(Mindex4)
colorbar
axis xy
yticks([0:15:60]);
yticklabels({'0','0.1','0.2','0.3','0.4'})
xticks([0:12:60]);
xticklabels({'0','0.2','0.4','0.6','0.8','1'})
hold off
grid on