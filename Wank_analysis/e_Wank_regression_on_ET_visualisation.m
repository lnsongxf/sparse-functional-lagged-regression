%{

This script continues in the analysis of the joint regression model and
produces the plots for the schematic figure on how the lagged regression
works (included in the paper).
The files "b_Wank_regression_on_E.m" and "c_Wank_regression_on_E.m" and
"d_Wank_regression_on_E.m" should have been run before.

Author: Tomas Rubin, tomas.rubin@gmail.com

Code developed on: MATLAB R2018a

%}



save_plots = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualisation of prediction


%% E withouth mean (substracted mu_est_ONB)

if save_plots, h_fig = figure('Position', [200,200,800, 150]); end

% without mean
adjustment = mu_est_ONB;
Eat_daily_centred = Eat_daily - mean(Eat_daily,1,'omitnan');

nDays = 7;
offset = 844; % 844
day_center = offset+nDays-3;
contribution_total = 0;
contributions = [];
for ii = 1:(nDays)
    t = offset+ii;
    
    if (t >= 1) && (t <= censor.nGridTime)
        
        subplot(1,7,ii)
        % substracted mean
        z = fSCBDegras( squeeze(krigingX_dyna_to_use.var(t,:,:)), onb.onbMatrix );
        lower_sim_dynamic = onb.onbMatrix*(krigingX_dyna_to_use.est(t,:)'-adjustment) - z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna_to_use.var(t,:,:))*onb.onbMatrix')));
        upper_sim_dynamic = onb.onbMatrix*(krigingX_dyna_to_use.est(t,:)'-adjustment) + z*sqrt(abs(diag( onb.onbMatrix*squeeze(krigingX_dyna_to_use.var(t,:,:))*onb.onbMatrix')));
        
        
        %         fill(24*[onb.gridSpace fliplr(onb.gridSpace)],[(lower_sim_dynamic') fliplr(upper_sim_dynamic')],...
        %             'b', 'FaceAlpha', 0.4) % dynamic band (simultaneous)
                
        
        % kriging visualisation
        plot( 24* onb.gridSpace, onb.onbMatrix * (krigingX_dyna.est( t,: )'-adjustment), "-k", "LineWidth", 1 )

        hold on
%         plot( 24* onb.gridSpace, upper_sim_dynamic, ":k", "LineWidth", 1 )
%         plot( 24* onb.gridSpace, lower_sim_dynamic, ":k", "LineWidth", 1 )
        % observed data
        nonzero_positions = find(Eat_indic(t,:));
        scatter( 24* onb.gridSpace(censor.onb.grid24(nonzero_positions)), Eat_daily_centred(t, nonzero_positions), 25, 'xk'  )
                
        
        
        hold off
        xlim([0 24])
        xticks([0 8 16 24])
        ylim([-150 200])
        
        if ii == 1
            ylabel({"observed atmospheric";"electricity E (V/m) (crosses);";...
                "estimated functional TS";" mean centred (solid line)"})
        end
    end
    
    % contribution from day "day_center"
    %     lag_real = -day_center + t;
    lag_real = +day_center - t;
    lag_index=lag_real + Tikh_transfer_1.numOfLags+1;
    if abs(lag_real)<= Tikh_transfer_1.numOfLags
        contribution = Tikh_transfer_1.ONB(lag_index,:) * ( krigingX_dyna_padded_to_use.est( t+dataFull_padding ,:)' - mu_est_ONB )    ;
        contributions = [contributions contribution];
        contribution_total = contribution_total + contribution;
        disp("t="+num2str(t)+", lag_real="+num2str(lag_real)+", lag_index="+num2str(lag_index)+", contribution="+num2str(contribution))
        
    end
    
    title("t="+num2str(t))
    
    %     title("$\hat{X}^{E}_{"+num2str(t)+"}$",'Interpreter','Latex')
end

contributions_E = contributions'

if save_plots
    set(h_fig,'Units','Inches');
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h_fig, 'figures_paper/Wank_E_848.pdf','-dpdf','-r0')
    close(h_fig)
end



%% time series T
contributions = [];
if save_plots, h_fig = figure('Position', [200,200,800,150]); end

for ii = 1:(nDays)
    subplot(1,7,ii)
    t = offset+ii;
    
    plot( 24*onb.gridSpace, onb.onbMatrix * (dataFullONB(t,:)' - T_mean), "-k", "LineWidth", 1)
    xlim([0 24])
	xticks([0 8 16 24])
    ylim([-15 15])
    
    % contribution from day "day_center"
    lag_real = +day_center - t;
    lag_index=lag_real + Tikh_transfer_2.numOfLags+1;
    if abs(lag_real)<= Tikh_transfer_2.numOfLags
        contribution = Tikh_transfer_2.ONB(lag_index,:) * ( krigingX_dyna_padded_to_use.est( t+dataFull_padding ,:)' - mu_est_ONB )    ;
        contributions = [contributions contribution];
        contribution_total = contribution_total + contribution;
        disp("t="+num2str(t)+", lag_real="+num2str(lag_real)+", lag_index="+num2str(lag_index)+", contribution="+num2str(contribution))
        
    end
    
    title("t="+num2str(t))
    if ii == 1, ylabel({"intraday temperature T (°C)";"mean centred";"yearly periodicity removed"}); end
end

contributions_T = contributions'

if save_plots
    set(h_fig,'Units','Inches');
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h_fig, 'figures_paper/Wank_T_848.pdf','-dpdf','-r0')
    close(h_fig)
end

%% filters E

if save_plots, h_fig = figure('Position', [200,200,800, 150]); end

for lag_real = -3:3
    subplot(1,7,-lag_real+4)
    plot( 24*onb.gridSpace, onb.onbMatrix * Tikh_transfer_1.ONB(lag_real + Tikh_transfer_1.numOfLags +1, :)', "-k" , "LineWidth", 1)
%     hline(0, ":k")
    ylim([-1.2 1.2])
    title("B^{(E)}_{"+num2str(lag_real)+"}" )
    xticks([0 8 16 24])    
    if lag_real == 3, ylabel({"estimated filter coefficients";"for atmoshperic electricity"}); end
end
pause(.1)

if save_plots
    set(h_fig,'Units','Inches');
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h_fig, 'figures_paper/Wank_filters_E.pdf','-dpdf','-r0')
    close(h_fig)
end


if save_plots, h_fig = figure('Position', [200,200,800, 150]); end

for lag_real = -3:3
    subplot(1,7,-lag_real+4)
    plot( 24*onb.gridSpace, onb.onbMatrix * Tikh_transfer_2.ONB(lag_real + Tikh_transfer_2.numOfLags +1, :)', "-k" , "LineWidth", 1)
%     hline(0, ":k")
    ylim([-4 4])
    title("B^{(\tau)}_{"+num2str(lag_real)+"}" )
    xticks([0 8 16 24])    
    if lag_real == 3, ylabel({"estimated filter coefficients";"for temperature"}); end
end
pause(.1)

if save_plots
    set(h_fig,'Units','Inches');
    pos = get(h_fig,'Position');
    set(h_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h_fig, 'figures_paper/Wank_filters_T.pdf','-dpdf','-r0')
    close(h_fig)
end


%% response
t = offset + (1:nDays)';
lag = fliplr(-3:3)';
response_that_day = response( t );
Tikh_ET_response_est_that_day = Tikh_ET_response_est( t )';
table( t, lag, contributions_E, contributions_T, response_that_day, Tikh_ET_response_est_that_day )
