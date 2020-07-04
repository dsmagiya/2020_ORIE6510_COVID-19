% SEIR Epidemic Model CTMC
% modified 4/10/2020
% updated 4/18/2020 Lucy Huo dh622
% updated 4/21/2020 dm852: add cumulative infected, denoted by III
% modified 5/7/2020 Desheng Ma dm852: add startday, endday of measures

%clear all 

function [sBenchmark,iBenchmark,ieBenchmark,eeBenchmark,rBenchmark,dBenchmark, iiiBenchmark]=SEIR_household_CTMC_Reopen()

TN=30000; % NYS total population is about 20*10^6, NYC total population is 8*10^6
eP=0.164; % elderly proportion
yP=0.20; % estimated youth/student proportion(<25) Cnesus gives up to <18 (~20%), but here would like to include university students


%first case in NYS 03/01/2020 
%Time = 45; % for CTMC only, calender time ->04/15/2020, death=214454
global From To T_measures

% set parameters
b = 0.08; % probability of transmission, 95% CI [0.05068, 0.05429]
k = [30 30 10]; % k = average # people exposed to by infected ppl
beta = k.*b; % infection rate, parameterize beta=k*b
incubation = 8;
delta = [1/inf 1/incubation 1/incubation]; % exposed -> infected rate, i.e. 1/incubation period
rec = 21; %three weeks to recover
gam = [1/rec 1/rec 1/rec]; % recover rate, 1/recovery period
alpha = [0.005 0.10 0.15]/5; % death rate

% initialize each group of people
s = {yP*TN (1-eP-yP)*TN-1 eP*TN}; % susceptible
i = {0 1 0}; % infected, initially 1 infected adult, NYS 1st confirmed case is adult
ie = {0 0 0}; % internally exposed
ee = {0 0 0}; % externally exposed
r = {0 0 0}; % recovered
d = {0 0 0}; % death
iii_ie = {0 0 0}; % cumulative infected from ie
iii_ee = {0 0 0}; % cumulative infected from ee

j = 1; % initialize time
t(1)=0; % initialize time
P = zeros; % initialize probability/intensities
N = TN; %initialize total population

sBenchmark=zeros(To,3);
iBenchmark=zeros(To,3);
ieBenchmark=zeros(To,3);
eeBenchmark=zeros(To,3);
rBenchmark=zeros(To,3);
dBenchmark=zeros(To,3);
iiiBenchmark = zeros(To,3); % 

%tic
while t(j)<To-From+1
    
    if t(j) > T_measures(1) && t(j) < T_measures(2)
        k = [30 30 10]/10; %k = average # people exposed to by infected ppl
        beta = k.*b; % infection rate, parameterize beta=k*b
    elseif t(j) > T_measures(2) && t(j) < T_measures(3)
        k = [30 30 10]; %k = average # people exposed to by infected ppl
        beta = k.*b; % infection rate, parameterize beta=k*b
    end
    
NN = N -(d{1}(j)+d{2}(j)+d{3}(j)); % deaths are removed from the population, sadly
    
% Intensities(rate)
P(1) = beta(1)*s{1}(j)*i{1}(j)/NN;
P(2) = beta(2)*s{2}(j)*i{2}(j)/NN;
P(3) = beta(3)*s{3}(j)*i{3}(j)/NN;
P(4) = (delta(1)/2)*ee{1}(j);
P(5) = (delta(2)/2)*ee{2}(j);
P(6) = (delta(3)/2)*ee{3}(j);
P(7) = (delta(1)/2)*ie{1}(j);
P(8) = (delta(2)/2)*ie{2}(j);
P(9) = (delta(3)/2)*ie{3}(j);
P(10) = (gam(1))*i{1}(j);
P(11) = (gam(2))*i{2}(j);
P(12) = (gam(3))*i{3}(j);
P(13) = alpha(1)*i{1}(j);
P(14) = alpha(2)*i{2}(j);
P(15) = alpha(3)*i{3}(j);
%P(16) = 1-P(1)-P(2)-P(3)-P(4)-P(5)-P(6)-P(7)-P(8)-P(9)-P(10)-P(11)-P(12)-P(13)-P(14)-P(15);
% Appeared in Chris' cide. Is this necessary? How to argue it?

% P/sum(P): normalize P to make entities less than 1
% dice rolling according to P
states = mnrnd(1,P/sum(P)); 
state = find(states);

% update entries for j+1 at the begining of each loop
       for k = 1:3
           s{k}(j+1) = s{k}(j);
           i{k}(j+1) = i{k}(j);
           ie{k}(j+1) = ie{k}(j);
           ee{k}(j+1) = ee{k}(j);
           r{k}(j+1) = r{k}(j);
           d{k}(j+1) = d{k}(j);
           iii_ee{k}(j+1) = iii_ee{k}(j);
           iii_ie{k}(j+1) = iii_ie{k}(j);
       end
       
% CTMC
% Gillespie algorithm
random = rand;
t(j+1)=-log(random)/sum(P)+t(j); % Time to next event using exponential distribution of lifetime

       if (state == 1) 
           % a kid is externally exposed, and the whole family is internally exposed immediately 
           nKids = binornd(5,0.3); 
           nAdult = 2; % a couple 
           nOld = binornd(5,0.3);
           
           if nKids > 0
           s{1}(j+1)=max(s{1}(j)-nKids,0);
           s{2}(j+1)=max(s{2}(j)-nAdult,0);
           s{3}(j+1)=max(s{3}(j)-nOld,0);
           
           ie{1}(j+1)=ie{1}(j)+min(nKids-1,s{1}(j));
           ie{2}(j+1)=ie{2}(j)+min(nAdult,s{2}(j));
           ie{3}(j+1)=ie{3}(j)+min(nOld,s{3}(j));
          
           ee{1}(j+1)=ee{1}(j)+1;
           end
           
       elseif (state == 2)
           % an adult is externally exposed, and whole family is internally exposed immediately 
           nKids = binornd(5,0.3);
           nAdult = 2; % a couple 
           nOld = binornd(5,0.3);
           
           if nAdult > 0
           s{1}(j+1) = max(s{1}(j)-nKids,0);
           s{2}(j+1) = max(s{2}(j)-nAdult,0);
           s{3}(j+1) = max(s{3}(j)-nOld,0);
           
           ie{1}(j+1) = ie{1}(j)+min(nKids,s{1}(j));
           ie{2}(j+1) = ie{2}(j)+min(nAdult-1,s{2}(j));
           ie{3}(j+1) = ie{3}(j)+min(nOld,s{3}(j));
           
           ee{2}(j+1) = ee{2}(j)+1;
           end
           
       elseif (state == 3) 
           % an old person is externally exposed, and whole family is internally exposed immediately 
           nKids = binornd(5,0.3);
           nAdult = 2; % a couple 
           nOld = binornd(5,0.3);
           
           if nOld > 0
           s{1}(j+1) = max(s{1}(j)-nKids,0);
           s{2}(j+1) = max(s{2}(j)-nAdult,0);
           s{3}(j+1) = max(s{3}(j)-nOld,0);
           
           ie{1}(j+1) = ie{1}(j)+min(nKids,s{1}(j));
           ie{2}(j+1) = ie{2}(j)+min(nAdult,s{2}(j));
           ie{3}(j+1) = ie{3}(j)+min(nOld-1,s{3}(j));
           
           ee{3}(j+1) = ee{3}(j)+1;
           end
           
       elseif (state == 4) 
           % ee kid is infected
           ee{1}(j+1) = ee{1}(j)-1; 
           i{1}(j+1) = i{1}(j)+1;
           
           iii_ee{1}(j+1) = iii_ee{1}(j)+1;
           
       elseif (state == 5) 
           % ee adult is infected
           ee{2}(j+1) = ee{2}(j)-1; 
           i{2}(j+1) = i{2}(j)+1; 
           
           iii_ee{2}(j+1) = iii_ee{2}(j)+1;
           
       elseif (state == 6) 
           % ee old is infected
           ee{3}(j+1) = ee{3}(j)-1; 
           i{3}(j+1) = i{3}(j)+1;
           
           iii_ee{3}(j+1) = iii_ee{3}(j)+1;
           
       elseif (state == 7) 
           % ie kid is infected
           ie{1}(j+1) = ie{1}(j)-1; 
           i{1}(j+1) = i{1}(j)+1;
           
           iii_ie{1}(j+1) = iii_ie{1}(j)+1;
           
       elseif (state == 8) 
           % ie adult is infected
           ie{2}(j+1) = ie{2}(j)-1; 
           i{2}(j+1) = i{2}(j)+1;
           
           iii_ie{2}(j+1) = iii_ie{2}(j)+1;
           
       elseif (state == 9) 
           % ie old is infected
           ie{3}(j+1) = ie{3}(j)-1; 
           i{3}(j+1) = i{3}(j)+1;
           
           iii_ie{3}(j+1) = iii_ie{3}(j)+1;
           
       elseif (state == 10) 
           % infected recovered
           i{1}(j+1) = i{1}(j)-1;
           r{1}(j+1) = r{1}(j)+1;
           
       elseif (state == 11) 
           % infected recovered
           i{2}(j+1) = i{2}(j)-1;
           r{2}(j+1) = r{2}(j)+1;
           
       elseif (state == 12) 
           % infected recovered
           i{3}(j+1) = i{3}(j)-1;
           r{3}(j+1) = r{3}(j)+1;
           
       elseif (state == 13)
           % infected died
           i{1}(j+1) = i{1}(j)-1;
           d{1}(j+1) = d{1}(j)+1;
           
       elseif (state == 14)
           % infected died
           i{2}(j+1) = i{2}(j)-1;
           d{2}(j+1) = d{2}(j)+1;
           
       elseif (state == 15)
           % infected died
           i{3}(j+1) = i{3}(j)-1;
           d{3}(j+1) = d{3}(j)+1;
          
       else
       end
       
       for kk = 1:To-From+1
           
           if t(j)>kk-1 && t(j)<=kk  
           sBenchmark(kk,:)=[s{1}(j+1),s{2}(j+1),s{3}(j+1)];
           iBenchmark(kk,:)=[i{1}(j+1),i{2}(j+1),i{3}(j+1)];
           ieBenchmark(kk,:)=[ie{1}(j+1),ie{2}(j+1),ie{3}(j+1)];
           eeBenchmark(kk,:)=[ee{1}(j+1),ee{2}(j+1),ee{3}(j+1)];
           rBenchmark(kk,:)=[r{1}(j+1),r{2}(j+1),r{3}(j+1)];
           dBenchmark(kk,:)=[d{1}(j+1),d{2}(j+1),d{3}(j+1)];
           iiiBenchmark(kk,:) = [iii_ee{1}(j+1)+iii_ie{1}(j+1),...
               iii_ee{2}(j+1)+iii_ie{2}(j+1),...
               iii_ee{3}(j+1)+iii_ie{3}(j+1)];
           end
           
       end
           
       j=j+1;
       
       
end %end for the while loop
%toc


for n = 2:To
    
    sBenchmark(n,sBenchmark(n,:)==0)=sBenchmark(n-1,sBenchmark(n,:)==0);
    iBenchmark(n,iBenchmark(n,:)==0)=iBenchmark(n-1,iBenchmark(n,:)==0);
    ieBenchmark(n,ieBenchmark(n,:)==0)=ieBenchmark(n-1,ieBenchmark(n,:)==0);
    eeBenchmark(n,eeBenchmark(n,:)==0)=eeBenchmark(n-1,eeBenchmark(n,:)==0);
    rBenchmark(n,rBenchmark(n,:)==0)=rBenchmark(n-1,rBenchmark(n,:)==0);
    dBenchmark(n,dBenchmark(n,:)==0)=dBenchmark(n-1,dBenchmark(n,:)==0);
    iiiBenchmark(n,iiiBenchmark(n,:)==0)=iiiBenchmark(n-1,iiiBenchmark(n,:)==0);
        
end


 %end for the function

 
%{
figure(1)

stairs(t,i{1}+i{2}+i{3},'b-','linewidth',1); 
hold on
%stairs(t,s{1}+s{2}+s{3},'m-','linewidth',1); 
stairs(t,ee{1}+ee{2}+ee{3},'g-','linewidth',1); 
stairs(t,ie{1}+ie{2}+ie{3},'y-','linewidth',1); 
stairs(t,r{1}+r{2}+r{3},'c-','linewidth',1); 
stairs(t,d{1}+d{2}+d{3},'r-','linewidth',1); 
xlim([0 Time]);
xlabel('Calender Time');
ylabel('Number of People');
%legend('Infected','Susceptible','Externally Exposed','Internally Exposed','Recovered','Death');
legend('Infected','Externally Exposed','Internally Exposed','Recovered','Death');
title('SEIR-household-CTMC');
grid on;

hold off
%}

%d{1}(end)+d{2}(end)+d{3}(end) % TotalDeath at End of Simulation
%i{1}(end)+i{2}(end)+i{3}(end) % TotalInfected at End of Simulation

end



