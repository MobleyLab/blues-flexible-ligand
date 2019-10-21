
NCMCstep = [300 600 1000 1500 2050 2500 3000];
nprop=5;
% calculating actual number of NCMC steps
NCMCactual = (0.4+0.6*nprop)*NCMCstep;
itrTotal=5000; % total number of iterations

for i = 1:length(NCMCstep)
    NCMCstep(i)
    % acceptance rate and error bar
    data = load( sprintf('acc_ncmc_%dNCMC.txt', NCMCstep(i)) );
    [ meanData(i), errData(i) ] = error(data);
end
meanData = meanData*100 %converting to percentages 
errData = errData*100 %converting to percentages

h=figure()
% plotting acceptance rate on left y-axis
yyaxis left
errorbar(NCMCactual, meanData, errData, 'o-', 'LineWidth', 1.2, 'MarkerFaceColor',[0,0.4470,0.7410]);
set(gca, 'FontSize',12);
xlabel('NCMC switching steps', 'FontSize',18);
ylabel('acceptance rate', 'FontSize',18);
ylim([0 25]);yticks(0:5:25);
hold on
% plotting moves accepted per million force evaluations on right y-axis
yyaxis right
moveData = meanData./(1000+NCMCactual)*10^6/100;
errMoveData = errData./(1000+NCMCactual)*10^6/100;
errorbar( NCMCactual, moveData , errMoveData, 's-', 'LineWidth', 1.2, 'MarkerFaceColor',[0.9,0.325,0.098]); 
ylabel('moves accepted/million f-ev', 'FontSize',18);
ylim([12 32]);
yticks(12:4:32)

ax = gca;
grid on
ax.GridColor = [0.1, 0.1, 0.1];


function [meanData, stdData] = error(data)
    % function for calculating mean and error bar for acceptance rates
    prob = data(:,1)./data(:,2);

    lenData = length(data);
    clear stdBlock

    for block = 2:20
        clear cumulativeStat probBlock
        lenBlock = floor(lenData/block);
        for i = 1:block
           cumulativeStat(i,:) = data(i*lenBlock, :);
           if i == 1
               probBlock(i) = cumulativeStat(i,1)/lenBlock;
           else
               probBlock(i) = (cumulativeStat(i,1)-cumulativeStat(i-1,1))/lenBlock;
           end
        end

        stdBlock(block-1) = std(probBlock)/sqrt(block);
    end
    %stdBlock
    stdData = max(stdBlock);
    meanData = data(end,1)/data(end,2);
end

