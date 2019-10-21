
NCMCstep = [300 600 1000 1500 2050 2500 3000];
nprop=5;
NCMCactual = (0.4+0.6*nprop)*NCMCstep;
itrTotal=5000;

for i = 1:length(NCMCstep)
    NCMCstep(i)
    % acceptance rate and error bar
    data = load( sprintf('../ligand1/NCMCstep_var/move-%dNCMC/acc_ncmc_%dNCMC.txt', NCMCstep(i), NCMCstep(i)) );
    [ meanData(i), errData(i) ] = error(data);
    
    % transition per million fev and error bar
    data = load( sprintf('../ligand1/NCMCstep_var/move-%dNCMC/torsion_lig1_MDNCMC.txt', NCMCstep(i)) );
    data(:,1) = data(:,1) * 0.001;
    
    iter = load( sprintf('../ligand1/NCMCstep_var/move-%dNCMC/iter-accp.txt', NCMCstep(i)) );
    fev = itrTotal*( 1000 + NCMCactual(i));
    iter = iter *0.002 ;
    [tr(i), errTr(i)] = transition(data, iter);
    tr(i) = tr(i)/fev*10^6;
    errTr(i) = errTr(i)/fev*10^6;
end
%converting acceptance rate to moves accepted / million f-ev
meanData
errData
meanData = meanData./(1000+NCMCactual)*10^6;
errData = errData./(1000+NCMCactual)*10^6;

h=figure()
yyaxis left
errorbar(NCMCactual, meanData, errData, 'o-', 'LineWidth', 1.2, 'MarkerFaceColor',[0,0.4470,0.7410]);
set(gca, 'FontSize',12);
ylabel('moves accepted/million f-ev', 'FontSize',18, 'position', [-1000, 21.5]);
ylim([8 32]);yticks(8:4:32);
hold on
yyaxis right
errorbar( NCMCactual, tr, errTr, 's-', 'LineWidth', 1.2, 'MarkerFaceColor',[0.9,0.325,0.098]); 


ax = gca;
grid on
ax.GridColor = [0.1, 0.1, 0.1];
xlabel('NCMC switching steps', 'FontSize',18, 'position', [6470, 5.5 ]);
ylabel('transitions/million f-ev', 'FontSize',18, 'position', [12800, 21.5]);
ylim([8 32]);
yticks(8:4:32)

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

