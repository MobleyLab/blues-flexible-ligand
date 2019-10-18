% this script is for plotting Figure 7 in paper. All the MD-NCMC simulations need to be performed ahead of time.
h = figure()

hold off;
x=0:10:10;
y=ones(1,length(x))*0.5;

% loading data
data = load('../ligand1/MD-NCMC-flip/torsion_lig1_MDNCMC.txt');
data(:,1) = data(:,1) * 0.001;
probLeft = getProb(data);
probRight = 1- probLeft;

% calculating probablity of state and associated error bar
time = data( length(data)/10: length(data)/10 : end , 1);
avgProbLeft = blockAvg ( probLeft, 10);
errorbar( time, avgProbLeft(:,1), avgProbLeft(:,2), 'o-', 'MarkerSize', 3, 'MarkerFaceColor',[0,0.4470,0.7410] );
hold on;
avgProbRight = blockAvg ( probRight, 10);

% plotting
errorbar( time, avgProbRight(:,1), avgProbRight(:,2), 's-','MarkerSize', 3, 'MarkerFaceColor',[0.9,0.325,0.098] );
plot(x,y,'-.','LineWidth', 1, 'Color',[0,0.4470,0.7410]);
ylim([-0.1 1.1]);
grid on;
xlim([0 10]);
legend( sprintf('orig   %0.2f',avgProbLeft(end,1) ),sprintf('flip    %0.2f',avgProbRight(end,1) ),'ideal 0.50');
set(gca, 'FontSize',14, 'YMinorTick', 'on', 'YMinorTick', 'on');

%%%calculationg moves accepted and transition per million fev
iter=load('../ligand1/MD-NCMC-flip/iter-accp.txt');
itrTotal = 5000;
NCMC=1000; nprop=5;
fev=itrTotal*( 1000 + 0.6*NCMC*nprop + 0.4*NCMC);
iter = iter *0.002 ; %converting to ns

[tr, errTr] = transition(data, iter);
tr/fev*10^6
errTr/fev*10^6

NCMCactual = 0.6*NCMC*nprop + 0.4*NCMC;
acc = length(iter)*10^6/itrTotal/(NCMCactual+1000);
title( {'{\bf3400 NCMC steps - flip move}', sprintf(' %0.0f moves accepted/million f-ev',acc), sprintf(' %0.1f transitions/million f-ev',tr/fev*10^6)}, 'FontWeight', 'Normal', 'position', [5.65, 1.1], 'FontSize', 14);
xlabel('time (in ns)', 'FontSize', 16);
ylabel('probability of state', 'FontSize', 16);
xticks(0:2.5:10); 


function [meanData, stdData] = error(data)
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

function [tr, errTr] = transition(data, iter)
    tr = 0;
    counter = []; %to keep track of where transitions are happening
    for i = 1:length(iter)
        idxData = round(iter(i) / 0.001) ;
        if data( idxData , 2) < -10
            if data( idxData + 1, 2) > -10
                tr = tr +1 ;
                counter(end+1) = idxData;
            end
        elseif data( idxData , 2) > -10
            if data( idxData + 1, 2) < -10
                tr = tr +1 ;
                counter(end+1) = idxData;
            end     
        end
    end
    
    %calculating error bar
    lenData = length(data);
    clear stdBlock

    for block = 2:20
        clear cumulativeStat probBlock
        lenBlock = floor(lenData/block);
        for i = 1:block
           a = ( counter >= lenBlock*(i-1));
           b = ( counter < lenBlock*i);
           combined = a + b;
           trBlock(i) = sum (combined == 2);
        end
        stdBlock(block-1) = std(trBlock)/sqrt(block);
    end
    %stdBlock
    errTr = max(stdBlock) ;
end
