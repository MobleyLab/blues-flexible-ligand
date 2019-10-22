# This script is for plotting fig 9 in the paper.

clear all
% from MD/NCMC, based on previous figures and analysis
leftProb(1) = 0.505;
errLeftProb(1) =0.005;

for i = 2:4
    if i == 2
      data = load('../ligand1/umbrellaSampling/umb_72win/PMF_umbrella.txt'); %72 windows
    elseif i == 3
      data = load('../ligand1/umbrellaSampling/umb_126win/PMF_umbrella.txt'); %126 windows
    else
      data = load('../ligand1/umbrellaSampling/umb_152win/PMF_umbrella.txt'); %152 widows
    end
    prob = exp(-data(:,2));
    prob = prob /sum(prob);
    data = data(prob>0.001,:); %only looking into values which have significant probability
    data(:,2)=exp(-data(:,2));
    data(:,3)=data(:,2).*data(:,3);

    left = data(data(:,1)<0,:);
    sumLeft = sum( left(:,2));
    errLeft = sqrt(sum(left(:,3).*left(:,3)));
    % right = data(data(:,1)>0,:);

    sumProb=sum(data(:,2));
    errProb=sqrt(sum(data(:,3).*data(:,3)));


    leftProb(i) = sumLeft/sumProb;
    errLeftProb(i) = leftProb(i) * sqrt( errLeft*errLeft/sumLeft/sumLeft + errProb*errProb/sumProb/sumProb);

end

h = figure()
total = [1,1,1,1];
h1 = bar(1:4, total, 0.6, 'FaceColor', [0.9, 0.9, 0.9]);
hold on;
h2 = bar(1:4, leftProb, 0.6, 'FaceColor', [0,0.4470,0.7410]);
errorbar(1:4, leftProb, errLeftProb,'.k', 'LineWidth', 1.5);

x=-1:5; y =0.5*ones(1,7);
plot(x,y,'k--','LineWidth', 1.2);
xlim([0.5 4.5]);
ylim([0 1.2]); yticks(0.2:0.2:1); set(gca, 'FontSize',16);

name={'MD/NCMC', 'Umb 72', 'Umb 126', 'Umb 152'};
set(gca,'xtick',[1:4],'xticklabel',name)
xtickangle(40);
ylabel('probability', 'FontSize',18, 'position', [0.1, 0.65]);

legend([h1(1) h2(1)],'flip','orig')
