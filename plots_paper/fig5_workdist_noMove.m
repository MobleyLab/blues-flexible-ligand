%%  gmx - alchemical, ring ocntaining R1, R2, R3-NVT %%%
%%%%%%%%%%%%
h=figure()
bins = -30:1:400;
acceptanceNVT = [0.0 15.8]; %for 1000 NCMC, obtained it from gmx.log file

% error = [0.0 1.5]; %for 1000 NCMC
for i = 1:2
    subplot(1,2,i);hold off;
    
    if i ==1
        data = load('../ligand1/MD-NCMC-noMove/wholeLig/work_ncmc_noMove_wholeLigand.txt');
    else 
        data = load('../ligand1/MD-NCMC-noMove/alchRegion/work_ncmc_noMove_alchRegion.txt');
    end
    i
    mean(data)
    std(data)
    if length(data) > 1000
       data = data(1:1000);
    end
    [x,y] = hist(-data, bins);
    x = x/length(data);
    bar(y,x,'FaceColor',[0,0.4470,0.7410]  );
    set(gca, 'FontSize',11);
    xlabel('\beta\itw', 'FontSize',15);
    ylabel('probability', 'FontSize',15);
    xlim([-15 100]);
    ylim([0 .13]);
    grid on;
end
