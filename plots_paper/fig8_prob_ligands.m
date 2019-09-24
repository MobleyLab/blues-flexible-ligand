h = figure()
% this script is for plotting Figure 8 in paper. All the MD-NCMC simulations need to be performed ahead of time.


folder={'ligand2', 'ligand3', 'ligand4'};
hold off;
x=0:10:100;
y=ones(1,length(x))*0.5;

for lig=1:2
    subplot(1,3,lig); hold off;
    data = load(sprintf('../ligand%d/MD-NCMC-flip/torsion_lig%d_MDNCMC.txt', i,i));
    data(:,1) = data(:,1) * 0.001;
    probLeft = getProb(data);
    probRight = 1- probLeft;
 
    noBlock = floor(length(data)/1000);   %1ns per time block
    time = data( length(data)/noBlock : length(data)/noBlock  : end , 1);
    avgProbLeft = blockAvg ( probLeft, noBlock );
    errorbar( time, avgProbLeft(:,1), avgProbLeft(:,2), 'o-', 'MarkerSize', 3, 'MarkerFaceColor',[0,0.4470,0.7410] );
    hold on;
    avgProbRight = blockAvg ( probRight, noBlock );
    errorbar( time, avgProbRight(:,1), avgProbRight(:,2), 's-','MarkerSize', 3, 'MarkerFaceColor',[0.9,0.325,0.098] );
    plot(x,y,'-.','LineWidth', 1, 'Color',[0,0.4470,0.7410]);
    ylim([-0.1 1.1]);
    grid on;
    xlim([0 noBlock]);
    legend( sprintf('orig   %0.2f',avgProbLeft(end,1) ),sprintf('flip    %0.2f',avgProbRight(end,1) ),'ideal 0.50');
    set(gca, 'FontSize',14, 'YMinorTick', 'on', 'YMinorTick', 'on');

    %%%calculationg moves accepted and transition per million fev
    if lig == 1
        iter=load('../ligand1/MD-NCMC-flip/iter-accp.txt');
        itrTotal = 5000;
    elseif lig == 2
        iter = load('../ligand2/MD-NCMC-flip/iter-accp.txt');
        itrTotal = 25000;
    end
    
    NCMC=1000;nprop=5;
    fev=itrTotal*( 1000 + 0.6*NCMC*nprop + 0.4*NCMC);
    iter = iter *0.002 ; %converting to ns
    [tr, errTr] = transition(data, iter)
    
    NCMCactual = 0.6*NCMC*nprop + 0.4*NCMC;
    acc = length(iter)*10^6/itrTotal/(NCMCactual+1000);
    tr/fev*10^6
    title( {sprintf('Ligand {\\bf%d}', (lig+1)), sprintf(' %0.1f \\pm %0.1f transitions/million f-ev',tr/fev*10^6, errTr/fev*10^6)}, 'FontWeight', 'Normal', 'FontSize', 14);
    xlabel('time (in ns)', 'FontSize', 16);
end
subplot(1,3,1)
ylabel('probability of state', 'FontSize', 18, 'position', [-1.5, 0.55]); xticks(0:2.5:10); 

subplot(1,3,2)
xticks(0:12.5:50);

% 
subplot(1,3,3)
hold off
data = load('../ligand4/MD-NCMC-flip/torsion_lig4_MDNCMC.txt');
data(:,1) = data(:,1) * 0.001;
probLeft = getProb(data);
probRight = 1- probLeft;

blocks = 20
time = data( length(data)/blocks: length(data)/blocks : end , 1);
avgProbRight = blockAvg ( probRight, blocks);
errorbar( time, avgProbRight(:,1), avgProbRight(:,2), 's-','LineWidth', 1, 'Color',[0,0.4470,0.7410], 'MarkerSize', 3, 'MarkerFaceColor',[0,0.4470,0.7410]);
hold on;
avgProbLeft = blockAvg ( probLeft, blocks);
errorbar( time, avgProbLeft(:,1), avgProbLeft(:,2), 'o-' ,'LineWidth', 1, 'Color',[0.9,0.325,0.098], 'MarkerSize', 3, 'MarkerFaceColor',[0.9,0.325,0.098] );

x=0:1:20;
y=ones(1,length(x))*0.6; 
plot(x,y,'-.','LineWidth', 1.5, 'Color',[0,0.4470,0.7410]);

y=ones(1,length(x))*0.4; 
hold on
plot(x,y,'-.','LineWidth', 1.5, 'Color',[0.9,0.325,0.098]);

iter = load('../ligand4/MD-NCMC-flip/iter-accp.txt');
itrTotal=20000;
NCMC=1000;nprop=5;
fev=itrTotal*( 1000 + 0.6*NCMC*nprop + 0.4*NCMC);
iter = iter *0.002 ; %converting to ns
[tr, errTr] = transition(data, iter)

NCMCactual = 0.6*NCMC*nprop + 0.4*NCMC;
acc = length(iter)*10^6/itrTotal/(NCMCactual+1000)
title( {sprintf('Ligand{\\bf%d}- MD/NCMC   ', 4), sprintf(' %0.1f \\pm %0.1f transitions/million f-ev',tr/fev*10^6,errTr/fev*10^6)}, 'FontWeight', 'Normal', 'FontSize', 14, 'position', [11,1.12]);

ylim([-0.1 1.1]);
grid on;
xlim([0 20]);
legend( sprintf('orig   %0.2f',avgProbRight(end,1) ),sprintf('flip    %0.2f',avgProbLeft(end,1) ),'orig exp 0.6','flip exp  0.4');
set(gca, 'FontSize',14, 'YMinorTick', 'on', 'YMinorTick', 'on');
xticks(0:5:20);
xlabel('time (in ns)', 'FontSize', 18, 'position', [10.5, -0.23]);

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

