%Muhammed Saadeddin Koçak 2232346
%% Clarification
clear
clc
close all
%% Parameter Definitions
propagationvelocity=0.5;%m/sec
startingpointofactivation=0;%starting point of activation shifted from the zsource in z direction
timebetweendepolarisationandrepolarisationwavefronts=0.3;%sec
strength=100;%mV
numberofsuccessivewavefronts=1;
timedelay=6;%sec
sigmai=0.4;
sigmao=1;
t=linspace(0,10,10001);%sec
sourcelength=5;%m
%% Locating Electrodes
E=[2.5 0.3 0.5;0.5 0.3 0.5;4.5 0.3 0.5];%z,r,theta
zsource=0;%m, z coordinate of initial point of source
d=zsource+startingpointofactivation;%location of depolarizing wave
r=zsource+startingpointofactivation;%location of repolarizing wave
%% Calculating Potentials
PD=zeros(length(t),3);
for i=1:length(t)
    for j=1:length(E)
        if d>=zsource+sourcelength
            PD(i,j)=0;
        else
            fun=@(x,y) 1./((sqrt((E(j,1)-d).^2+(E(j,2)-x).^2+(E(j,3)-y).^2)).^3);
            polarfun=@(theta,r) fun(r.*cos(theta),r.*sin(theta));
            q=integral2(polarfun,0,2*pi,0,0.05);
            PD(i,j)=q*sigmai*strength*(E(j,1)-d)/(4*pi*sigmao);
        end
    end
    d=d+propagationvelocity*0.001;
end
PR=zeros(length(t),3);
for i=timebetweendepolarisationandrepolarisationwavefronts*1000:length(t)
    for j=1:length(E)
        if r>=zsource+sourcelength
            PR(i,j)=0;
        else
            fun=@(x,y) 1./((sqrt((E(j,1)-r).^2+(E(j,2)-x).^2+(E(j,3)-y).^2)).^3);
            polarfun=@(theta,r) fun(r.*cos(theta),r.*sin(theta));
            q=integral2(polarfun,0,2*pi,0,0.05);
            PR(i,j)=-q*sigmai*strength*(E(j,1)-r)/(4*pi*sigmao);
        end
    end
    r=r+propagationvelocity*0.001;
end
if numberofsuccessivewavefronts==2
    d=zsource+startingpointofactivation;%location of depolarizing wave
    r=zsource+startingpointofactivation;%location of repolarizing wave
    PD2=zeros(length(t),3);
    for i=timedelay*1000:length(t)
        for j=1:length(E)
            if d>=zsource+sourcelength
                PD2(i,j)=0;
            else
                fun=@(x,y) 1./((sqrt((E(j,1)-d).^2+(E(j,2)-x).^2+(E(j,3)-y).^2)).^3);
                polarfun=@(theta,r) fun(r.*cos(theta),r.*sin(theta));
                q=integral2(polarfun,0,2*pi,0,0.05);
                PD2(i,j)=q*sigmai*strength*(E(j,1)-d)/(4*pi*sigmao);
            end
        end
        d=d+propagationvelocity*0.001;
    end
    PD=PD+PD2;
    PR2=zeros(length(t),3);
    for i=(timebetweendepolarisationandrepolarisationwavefronts+timedelay)*1000:length(t)
        for j=1:length(E)
            if r>=zsource+sourcelength
                PR2(i,j)=0;
            else
                fun=@(x,y) 1./((sqrt((E(j,1)-r).^2+(E(j,2)-x).^2+(E(j,3)-y).^2)).^3);
                polarfun=@(theta,r) fun(r.*cos(theta),r.*sin(theta));
                q=integral2(polarfun,0,2*pi,0,0.05);
                PR2(i,j)=-q*sigmai*strength*(E(j,1)-r)/(4*pi*sigmao);
            end
        end
        r=r+propagationvelocity*0.001;
    end
    PR=PR+PR2;
end
P=PD+PR;
V1=P(:,1)-P(:,2);
V2=P(:,2)-P(:,3);
V3=P(:,3)-P(:,1);
%% Presenting Results
figure
hold on
plot(t,V1(:,1))
plot(t,V2(:,1))
plot(t,V3(:,1))
xlabel('Time(sec)')
ylabel('Potential Difference(mV)')
legend('P1-P2','P2-P3','P3-P1')

figure
hold on
plot(t,P(:,1))
plot(t,P(:,2))
plot(t,P(:,3))
xlabel('Time(sec)')
ylabel('Potential(mV)')
legend('Electrode 1','Electrode 2','Electrode 3')

figure
hold on
plot(t,PD(:,1))
plot(t,PR(:,1))
plot(t,P(:,1))
xlabel('Time(sec)')
ylabel('Potential(mV)')
legend('Depolarizing Pulse','Repolarizing Pulse','Potential Measured')
% legend('Electrode 1','Electrode 2','Electrode 3')