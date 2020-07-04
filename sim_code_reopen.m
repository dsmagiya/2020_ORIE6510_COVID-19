% SEIR Epidemic Model CTMC
% modified 4/10/2020
% updated 4/18/2020 Lucy Huo dh622
% updated 4/21/2020 Desheng Ma dm852: add cumulative infected, denoted by III
% modified 5/7/2020 Desheng Ma dm852: add startday, endday of measures

clear all;
tic

figure(2)

startday = 5:5:20;
% endday = 50:10:80;
% startday = 10;
endday = 70;
% for day = endday
for seir = 1:4
    
f=waitbar(0,'SEIR');

global From To T_measures

From = 1; % stupid mistake here, From has to be always 1
To = 100;
T_measures = [startday(seir) endday 100];

simRound=100;
sCum=zeros(simRound,To-From+1,3);
iCum=zeros(simRound,To-From+1,3);
ieCum=zeros(simRound,To-From+1,3);
eeCum=zeros(simRound,To-From+1,3);
rCum=zeros(simRound,To-From+1,3);
dCum=zeros(simRound,To-From+1,3);
iiiCum = zeros(simRound,To-From+1,3);

for n=1:simRound
    [s,i,ie,ee,r,d,iii]=SEIR_household_CTMC_Reopen();
    sCum(n,:,:)=s;
    iCum(n,:,:)=i;
    ieCum(n,:,:)=ie;
    eeCum(n,:,:)=ee;
    rCum(n,:,:)=r;
    dCum(n,:,:)=d;
    iiiCum(n,:,:)=iii;
    waitbar(n/simRound,f,'SEIR');
end
close(f);
avgDailyS=[20*10^8;mean(sum(sCum,3))'];
avgDailyI=[1;mean(sum(iCum,3))'];
avgDailyIE=[0;mean(sum(ieCum,3))'];
avgDailyEE=[0;mean(sum(eeCum,3))'];
avgDailyR=[0;mean(sum(rCum,3))'];
avgDailyD=[0;mean(sum(dCum,3))'];
avgDailyIII=[0;mean(sum(iiiCum,3))'];

stdDailyS=std(sum(sCum,3))';
stdDailyI=std(sum(iCum,3))';
stdDailyIE=std(sum(ieCum,3))';
stdDailyEE=std(sum(eeCum,3))';
stdDailyR=std(sum(rCum,3))';
stdDailyD=std(sum(dCum,3))';
stdDailyIII=std(sum(iiiCum,3))';

z=norminv(0.975);
avgDailySCI=[20*10^8,20*10^8;avgDailyS(2:end)-(z/sqrt(simRound)).*stdDailyS,avgDailyS(2:end)+(z/sqrt(simRound)).*stdDailyS];
avgDailyICI=[1,1;avgDailyI(2:end)-(z/sqrt(simRound)).*stdDailyI,avgDailyI(2:end)+(z/sqrt(simRound)).*stdDailyI];
avgDailyIECI=[0,0;avgDailyIE(2:end)-(z/sqrt(simRound)).*stdDailyIE,avgDailyIE(2:end)+(z/sqrt(simRound)).*stdDailyIE];
avgDailyEECI=[0,0;avgDailyEE(2:end)-(z/sqrt(simRound)).*stdDailyEE,avgDailyEE(2:end)+(z/sqrt(simRound)).*stdDailyEE];
avgDailyRCI=[0,0;avgDailyR(2:end)-(z/sqrt(simRound)).*stdDailyR,avgDailyR(2:end)+(z/sqrt(simRound)).*stdDailyR];
avgDailyDCI=[0,0;avgDailyD(2:end)-(z/sqrt(simRound)).*stdDailyD,avgDailyD(2:end)+(z/sqrt(simRound)).*stdDailyD];
avgDailyIIICI=[0,0;avgDailyIII(2:end)-(z/sqrt(simRound)).*stdDailyIII,avgDailyIII(2:end)+(z/sqrt(simRound)).*stdDailyIII];

%%%%%%%%%%%%%%%%%%%%%%%
% load('cleanNYSData.mat')

%%% Simulated SEIR-household-CTMC
subplot(2,2,seir)

stairs((From-1:To)',avgDailyI,'m-','linewidth',1); 
hold on
% stairs((From-1:To)',avgDailyS,'b-','linewidth',1);
stairs((From-1:To)',avgDailyIE,'g-','linewidth',1); 
stairs((From-1:To)',avgDailyEE,'b-','linewidth',1); 
stairs((From-1:To)',avgDailyR,'c-','linewidth',1); 
stairs((From-1:To)',avgDailyD,'r-','linewidth',1); 
stairs((From-1:To)',avgDailyIII,'k-','linewidth',1); 

xlim([From-1 To]);
ylim([0 30000]);
xlabel('Calender Time');
ylabel('Number of People');
%legend('Infected','Susceptible','Internally Exposed','Externally Exposed','Recovered','Death');
legend('Infected','Internally Exposed','Externally Exposed','Recovered','Death','Cumulative Infected');
stitle = sprintf('SEIR-household-CTMC measures: day %i to %i', startday(seir), endday);
title(stitle)
grid on;

hold off

end

toc
% 
% %%% SEIR-household-CTMC Recovered Comparison
% figure(2)
% stairs((From-1:To)',avgDailyR,'c-','linewidth',1); 
% hold on
% stairs((From-1:To)',avgDailyRCI(:,1),'r-','linewidth',1); 
% stairs((From-1:To)',avgDailyRCI(:,2),'r-','linewidth',1); 
% stairs((From-1:To)',[0;cleanNYSData.Recovered(From:To)],'y-','linewidth',1); 
% 
% xlim([From-1 To]);
% xlabel('Calender Time');
% ylabel('Number of People');
% legend('Average Recovered','Upper 95% CI','Lower 95% CI','Actual Recovered');
% title('SEIR-household-CTMC Recovered Comparison');
% grid on;
% 
% hold off
% 
% %%% SEIR-household-CTMC Death Comparison
% figure(3)
% stairs((From-1:To)',avgDailyD,'c-','linewidth',1); 
% hold on
% stairs((From-1:To)',avgDailyDCI(:,1),'r-','linewidth',1); 
% stairs((From-1:To)',avgDailyDCI(:,2),'r-','linewidth',1); 
% stairs((From-1:To)',[0;cleanNYSData.Deaths(From:To)],'y-','linewidth',1); 
% 
% xlim([From-1 To]);
% xlabel('Calender Time');
% ylabel('Number of People');
% legend('Average Deaths','Upper 95% CI','Lower 95% CI','Actual Deaths');
% title('SEIR-household-CTMC Deaths Comparison');
% grid on;
% 
% 
% %%% SEIR-household-CTMC Confirmed(/Infected) Comparison
% figure(4)
% stairs((From-1:To)',avgDailyIII,'c-','linewidth',1); 
% hold on
% stairs((From-1:To)',avgDailyIIICI(:,1),'r-','linewidth',1); 
% stairs((From-1:To)',avgDailyIIICI(:,2),'r-','linewidth',1); 
% stairs((From-1:To)',[0;cleanNYSData.Confirmed(From:To)],'y-','linewidth',1); 
% 
% xlim([From-1 To]);
% xlabel('Calender Time');
% ylabel('Number of People');
% legend('Average Infected','Upper 95% CI','Lower 95% CI','Actual Confirmed');
% title('SEIR-household-CTMC Confirmed(/Infected) Comparison');
% % grid on;
% hold off
