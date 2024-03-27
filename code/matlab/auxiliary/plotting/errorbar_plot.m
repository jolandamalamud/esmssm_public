function errorbar_plot(data, xlabelling, ylabelling, newfig)

    if newfig
        
        figure('DefaultAxesFontSize',18);
    end
        
        bar(nanmean(data), 'grouped');
        hold on;
        er = errorbar(nanmean(data), nanstd(data));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        ylabel(ylabelling);
        set(gca, 'FontSize', 18, 'xticklabel', xlabelling, 'xticklabelRotation', 45);
       

end