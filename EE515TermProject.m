%% Muhammed Saadeddin Koçak
%% Clear
clear
%% Initialization
%Remember to update variables according to project description
%Need further variable definitions
a=0.03;%axon radius,cm
Ri=94;%intracellular resistivity,ohmcm
ri=Ri/(pi*a^2);%intracellular resistance per unit length,ohm/cm
Re=20;%extracellular resistivity,ohmcm
re=Re/(pi*a^2);%extracellular resistance per unit length,ohm/cm
Cm=1;%uF/cm^2
EK=-102;%mV
ENa=25;%mV
EL=-79.387;%mV
gKmax=36;%mS/cm^2
gNamax=120;%mS/cm^2
gLmax=0.3;%mS/cm^2
t=linspace(-99,15000,15100);%usec
l=linspace(0,1,101);%cm
n=zeros(length(t),length(l));
m=zeros(length(t),length(l));
h=zeros(length(t),length(l));
an=zeros(length(t),length(l));
bn=zeros(length(t),length(l));
am=zeros(length(t),length(l));
bm=zeros(length(t),length(l));
ah=zeros(length(t),length(l));
bh=zeros(length(t),length(l));
Vm=zeros(length(t),length(l));
gK=zeros(length(t),length(l));
gNa=zeros(length(t),length(l));
Is=zeros(length(t),length(l));
% gK(:,1)=gKmax*(n(:,1).^4);
% gNa(:,1)=gNamax*(m(:,1).^3).*(h(:,1));
% Vr=(EL*gLmax+EK*gNa(1,1)+EK*gK(1,1))/(gLmax+gNa(1,1)+gK(1,1));%See pg. 176
Vr=-90;%mV
Vm(1:100,:)=Vr;%mV
vm=zeros(length(t),length(l));
Im=zeros(length(t),length(l));
impulseduration=200;%usec
numberofconsecutivestimuli=2;
timedelaybetweenstimuli=8000;%usec
impulseamplitude=500;%uA/cm^2
%% Loops
figure;
for i=100:length(t)
    for j=1:length(l)
           if t(i)>-1 && t(i)<=impulseduration
               Is(i,2)=impulseamplitude;%uA/cm^2
           else
               Is(i,2)=0;
           end
           %For the second stimuli
           if numberofconsecutivestimuli~=1
           if t(i)>timedelaybetweenstimuli && t(i)<=timedelaybetweenstimuli+impulseduration
               Is(i,2)=impulseamplitude;%uA/cm^2
           end
           end
%         Im(i,j)=(((Vm(i,j-1)-2*Vm(i,j)+Vm(i,j+1))-re*Is(i,j))/(2*pi*a*(ri+re)))+Is(i,j);%According
        %to book description
        if j~=1 && j~=length(l)
        Im(i,j)=((Vm(i,j-1)-2*Vm(i,j)+Vm(i,j+1))*a/((2*Ri)*(0.01^2)))+Is(i,j);
        end
        %Calculating ionic currents
        an(i,j)=((0.01*(10-vm(i,j)))/(exp((10-vm(i,j))/10)-1))/1000;%1/usec
        bn(i,j)=(0.125*exp((-vm(i,j))/80))/1000;%1/usec
        am(i,j)=((0.1*(25-vm(i,j)))/(exp(0.1*(25-vm(i,j)))-1))/1000;%1/usec
        bm(i,j)=(4*exp((-vm(i,j))/18))/1000;%1/usec
        ah(i,j)=(0.07*exp((-vm(i,j))/20))/1000;%1/usec
        bh(i,j)=(1/(exp((30-vm(i,j))/10)+1))/1000;%1/usec
        n(100,j)=an(100,j)/(an(100,j)+bn(100,j));
        m(100,j)=am(100,j)/(am(100,j)+bm(100,j));
        h(100,j)=ah(100,j)/(ah(100,j)+bh(100,j));
        gK(i,j)=gKmax*(n(i,j).^4);%mS/cm^2
        gNa(i,j)=gNamax*(m(i,j).^3)*h(i,j);%mS/cm^2
        IK(i,j)=gK(i,j)*(Vm(i,j)-EK);%uA/cm^2
        INa(i,j)=gNa(i,j)*(Vm(i,j)-ENa);%uA/cm^2
        IL(i,j)=gLmax*(Vm(i,j)-EL);%uA/cm^2
        deltan(i,j)=(an(i,j)*(1-n(i,j)))-(bn(i,j)*n(i,j));
        deltam(i,j)=(am(i,j)*(1-m(i,j)))-(bm(i,j)*m(i,j));
        deltah(i,j)=(ah(i,j)*(1-h(i,j)))-(bh(i,j)*h(i,j));
        deltaVm(i,j)=(Im(i,j)-IK(i,j)-INa(i,j)-IL(i,j))*10/(Cm*1000);%mV
        n(i+1,j)=n(i,j)+deltan(i,j);
        m(i+1,j)=m(i,j)+deltam(i,j);
        h(i+1,j)=h(i,j)+deltah(i,j);
        Vm(i+1,j)=Vm(i,j)+deltaVm(i,j);
        vm(i+1,j)=vm(i,j)+deltaVm(i,j);
    end
    plot(l,Vm(i,1:length(l)));
    title('Membrane Voltage vs Position');
    xlabel('Position (cm)');
    ylabel('Voltage(mV)');
    drawnow;
end
%% Results
% plot(t,Vm(1:length(t),21));
% title('Membrane Voltage vs Time at 0.2cm');
% xlabel('Time (usec)');
% ylabel('Voltage(mV)');
figure
hold on
plot(t,Vm(1:length(t),21));
plot(t,Vm(1:length(t),41));
plot(t,Vm(1:length(t),61));
plot(t,Vm(1:length(t),81));
legend({'0.2cm','0.4cm','0.6cm','0.8cm'});
title('Membrane Voltage vs Time at Various Positions');
xlabel('Time (10usec)');
ylabel('Voltage(mV)');
hold off


% figure
% plot(l,Vm(1000,1:length(l)));
% title('Position vs Membrane Voltage at t=10000us');
% xlabel('Position (cm)');
% ylabel('Voltage(mV)');
figure
hold on
plot(l,Vm(3000,1:length(l)));
plot(l,Vm(6000,1:length(l)));
plot(l,Vm(9000,1:length(l)));
plot(l,Vm(12000,1:length(l)));
plot(l,Vm(15000,1:length(l)));
legend({'30000usec','60000usec','90000usec','120000usec','150000usec'});
title('Membrane Voltage vs Position at Various Time Instants');
xlabel('Position (cm)');
ylabel('Voltage(mV)');
hold off

figure
hold on
plot(l,h(3000,1:length(l)));
plot(l,h(9000,1:length(l)));
plot(l,h(15000,1:length(l)));
plot(l,n(3000,1:length(l)));
plot(l,n(9000,1:length(l)));
plot(l,n(15000,1:length(l)));
plot(l,m(3000,1:length(l)));
plot(l,m(9000,1:length(l)));
plot(l,m(15000,1:length(l)));

legend({'h at 30000usec','h at 90000usec','h at 150000usec','n at 30000usec','n at 90000usec','n at 150000usec','m at 30000usec','m at 90000usec','m at 150000usec'});
title('Activation Parameters vs Position at Various Time Instants');
xlabel('Position (cm)');
ylabel('Voltage(mV)');
hold off