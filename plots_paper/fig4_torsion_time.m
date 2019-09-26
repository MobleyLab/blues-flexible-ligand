% h=figure()
hold off
%%% plotting MD
h1=subplot(2,3,1:2)
%written after every 10 fs
tr = 1;
fev=100*1000/0.002;
data = load('../ligand1/MD/torsion_lig1_MD.txt');
data(:,1) = data(:,1) * 0.01; 
plot(data(:,1), data(:,2),'.')
%plot(data(1:2000,1), data(4001:6000,2),'.')
grid on; set(gca, 'FontSize',11);
ylabel('torsion angle', 'FontSize',16);
xlabel('time (in ns)', 'FontSize',16);
title(sprintf('MD, %0.2f transitions/million f-ev', tr/fev*10^6), 'position',[55 205], 'FontSize',14);



%%%plotting MD probability of convergence
h2=subplot(2,3,3); hold off;
x=0:10:110;
y=ones(1,length(x))*0.5;
probLeft = getProb(data);
probRight = 1- probLeft;

time = data( length(data)/10: length(data)/10 : end , 1);
avgProbLeft = blockAvg ( probLeft, 10);
errorbar( time, avgProbLeft(:,1), avgProbLeft(:,2), 'o-' ,'LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor',[0,0.4470,0.7410] );
hold on;
avgProbRight = blockAvg ( probRight, 10);
errorbar( time, avgProbRight(:,1), avgProbRight(:,2), 's-','LineWidth', 1, 'MarkerSize', 3, 'MarkerFaceColor',[0.9,0.325,0.098]);
plot(x,y,'-.','LineWidth', 1, 'Color',[0,0.4470,0.7410]);
ylim([-0.1 1.1]);
grid on;
xlim([0 100]); set(gca, 'FontSize',11);
legend( sprintf('orig   %0.2f',avgProbLeft(end,1) ),sprintf('flip    %0.2f',avgProbRight(end,1) ),'ideal 0.50');
title( 'MD', 'FontSize',14 );
xlabel('time (in ns)', 'FontSize',16); ylabel('probability of state', 'FontSize',16);
xticks(0:25:100); 

% plotting NCMC, 5000 itr, 1000NCMC, 1000MD
h3=subplot(2,3,4:5); hold off;
data = load('../ligand1/MD-NCMC/torsion_lig1_MDNCMC.txt');
data(:,1) = data(:,1) * 0.001;
plot(data(:,1), data(:,2),'.')
hold on

iter=load('../ligand1/MD-NCMC/iter-accp.txt');

nprop=5; NCMC=1000; itrTotal=5000;
fev=itrTotal*( 1000 + 0.6*NCMC*nprop + 0.4*NCMC);
iter = iter *0.002 ; %converting to ns
y=-200:20:200;
tr = 0;
for i = 1:length(iter)
    idxData = round(iter(i) / 0.001) ;
    if data( idxData , 2) < -10
        if data( idxData + 1, 2) > -10
            tr = tr +1 ;
        end
    elseif data( idxData , 2) > -10
        if data( idxData + 1, 2) < -10
            tr = tr +1 ;
        end     
    end
end
tr
set(gca, 'FontSize',11);
ylabel('torsion angle', 'FontSize',16);
xlabel('time (in ns)', 'FontSize',16); 
title(sprintf('MD/NCMC, %0.2f transitions/million f-ev', tr/fev*10^6 ), 'position',[5.5 205], 'FontSize',14);
grid on;

%%% plot MD/NCMC probability of convergence
h4=subplot(2,3,6); hold off;
x=0:10:100;
y=ones(1,length(x))*0.5;
probLeft = getProb(data);
probRight = 1- probLeft;

time = data( length(data)/10: length(data)/10 : end , 1);
avgProbLeft = blockAvg ( probLeft, 10);
errorbar( time, avgProbLeft(:,1), avgProbLeft(:,2), 'o-', 'LineWidth', 1,'MarkerSize', 3, 'MarkerFaceColor',[0,0.4470,0.7410 ]);
hold on;
avgProbRight = blockAvg ( probRight, 10);
errorbar( time, avgProbRight(:,1), avgProbRight(:,2), 's-','LineWidth', 1,'MarkerSize', 3, 'MarkerFaceColor',[0.9,0.325,0.098] );
plot(x,y,'-.','LineWidth', 1, 'Color',[0,0.4470,0.7410]);
ylim([-0.1 1.1]);
grid on;
xlim([0 10]);

totaltime= time(end);
set(gca, 'FontSize',11); title('MD/NCMC' , 'FontSize',14); 
legend( sprintf('orig   %0.2f',avgProbLeft(end,1) ),sprintf('flip    %0.2f',avgProbRight(end,1) ),'ideal 0.50');
xticks(0:2.5:10); xlabel('time (in ns)', 'FontSize',16); ylabel('probability of state', 'FontSize',16); 


