% clear all 
%% function used to plot the figures found in the article "An infrastructure perspective for\\
%% enhancing multi-functionality of forests:A conceptual modeling approach" submitted to the Earths' future journal
%% The function takes inputs myplot(construction investment ratio, tax imposition parameter, fixed initial value for bigtrees,\\
%% fixed initial value for small trees, value for tree removal, value for deadwood removal, a time horizon for the model simlation) Figures (2-5)
function myplot(CONSINV,TAX,STRATA1,STRATA2,kk,oo,TTT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% Simulation of points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ooo=[];
iii=[];
NN=[];
PP=[];
AM=[];
XX=1;
AC=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Model simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=CONSINV;
dt=0.1;

%% simulation for a model evolution for time
for alpha1=ee

EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
TCT=TAX*0.5;
TCF=TAX;
EE=[];
PP=[];
CC=[];
TT=[];
MM=[];
MMM=[];
CCT=[];
CCF=[];
wood=[];
BB1=[];
BB2=[];
minstrata1=[];
testt=[];
MFII=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y1(1)=Q*STRATA1;y2(1)=Q*STRATA2;y3(1)=25*Q;y4(1)=0.3;y5(1)=0.3;y6(1)=100;y7(1)=0.1;BA(1)=0.07;E(1)=0.01;P(1)=0.001;II=0.5;hh=kk*Q;o=oo*Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                 for t=2:TTT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% forest dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Strata 1 - big trees
                        if y1(t-1)>0
                            hh=kk*Q;
                            y1(t) = y1(t-1) + dt*(h*y2(t-1)*(1-u*g1*y1(t-1)) - y1(t-1)*d - hh*HH(y5(t-1),y4(t-1),SI0,SIm)/(v1));
                        else
                            y1(t)=0 ;
                            hh=0;
                            II=0;
                        end 
                %%% Strata 2 - small trees
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
                            o=0;  
                        end
                %%% Annual budget
                        BA(t)= dt*(U*(TCF*(hh*(HH(y5(t-1),y4(t-1),SI0,SIm))*(p-cm)+...
                                   (o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))*(pde-cde))+gamma+TCT*y6(t-1)));  
                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Infrastructures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Infrastructure construction dynamics
                        if y4(t-1)>0
                            y4(t) = y4(t-1) +   dt*(II*(alpha1)*BA(t)*u1)*mu*y4(t-1)*(1-(y4(t-1)/1)); 
                        else
                            y4(t)=0;
                        end 

                %%% Infrastructures state dynamics
                            ez=((1-alpha1)*II*BA(t)*u2*(1/(y4(t-1)+k)));
                            y5(t) = y5(t-1) + dt*((ez)*y5(t-1)*((1-y5(t-1))/1)+(y4(t)-y4(t-1))*((SIm-y5(t-1))/y4(t-1))- delta*y5(t-1)) ; %%logistic equation
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
                                alphaC=6;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            else
                                alphaC=6;
                                y6(t) = y6(t-1) + dt*(y6(t-1)*(E + HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT) - alphaT*y6(t-1)- alphaC));
                            end
                        else
                            y6(t)=0;
                        end
                    
                    
                    
%                     EE=[EE,E];
%                     PP=[PP,HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT)];
%                     CCF=[CCF,HH(y5(t-1),y4(t-1),SI0,SIm)];
%                     CC=[CC,-alphaT*y6(t-1)];
%                     TT=[TT,t];
%                     MM=[MM,II*(alpha1*BA(t)*u1*(1/(y4(t-1)+k)))];
%                     MMM=[MMM,-delta*y5(t-1)];
%                     CCT=[CCT,(aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2))];
%                     wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
%                     BB1=[BB1,dt*(TCF*(hh*(HH(y5(t-1),y4(t-1),SI0,SIm))*(p-cm)+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))*(pde-cde)))];
%                     BB2=[BB2,dt*TCT*y6(t-1)];
%                     testt=[testt,aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2)];
%                  end
%                     EE=[EE,E];
%                     PP=[PP,HHT(y4(t-1),y5(t-1),CI0,SIm1,SI01,aT,CT)];
%                     CCF=[CCF,HH(y5(t-1),y4(t-1),SI0,SIm)];
%                     CC=[CC,-alphaT*y6(t-1)];
%                     TT=[TT,t];
%                     MM=[MM,II*(alpha1*BA(t)*u1*(1/(y4(t-1)+k)))];
%                     MMM=[MMM,-delta*y5(t-1)];
%                     CCT=[CCT,(aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2))];
%                     wood=[wood,hh*(HH(y5(t-1),y4(t-1),SI0,SIm))+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))];
%                     BB1=[BB1,dt*(TCF*(hh*(HH(y5(t-1),y4(t-1),SI0,SIm))*(p-cm)+(o*v1*y1(t-1)*d*pe*pa*HH(y5(t-1),y4(t-1),SI0,SIm))*(pde-cde)))];
%                     BB2=[BB2,dt*TCT*y6(t-1)];
%                     testt=[testt,aT*exp((-1/2)*((y4(t-1)-CI0)/CT)^2)];
%                     MFII=[MFII,0.33*(y6(:)/max(y6(:)))+0.33*(wood(:)/max(wood(:)))+0.33*(y3(:)/max(y3(:)))];
end







%% not necessary
% set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
figure(12222221)
hold on
grid on
subplot(2,3,1)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(y1,'lineWidth',2)
xlabel('Time in years')
ylabel('Number of large trees/ha')
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,2)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(TT,wood,'lineWidth',2)
xlabel('Time in years')
ylabel('Wood removal volume (m^3/ha)')
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,3)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(y6,'lineWidth',2)
xlabel('Time in years')
ylabel('Number of tourists/ha')
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,4)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(y3,'lineWidth',2)
xlabel('Time in years')
ylabel('Deadwood volume (m^3/ha)')
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,5)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(y5,'lineWidth',2)
xlabel('Time in years')
ylabel('State of roads (S_I)')
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,6)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(y4,'lineWidth',2)
xlabel('Time in years')
ylabel('Quantity of roads (C_I)')
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
grid on
hold off



figure(111112121)
hold on
grid on
subplot(2,3,1)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(BA,'lineWidth',2)
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,2)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(TT,BB1,'lineWidth',2)
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
hold off
hold on
grid on
subplot(2,3,3)
set(gca,'linestyleorder',{'-','--','-.'}, 'colororder',[0 0 0],'nextplot','add')
plot(TT,BB2,'lineWidth',2)
xt = get(gca, 'XTick');    
set(gca, 'XTick', xt, 'XTickLabel', xt/10) 
grid on
hold off






end


