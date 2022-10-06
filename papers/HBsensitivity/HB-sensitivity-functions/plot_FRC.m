function plot_FRC(frc_data, items2plot, plotOptions)
%DESCRIPTION:
% void function to create personalized plots of frcs
% 
% INPUTS:
% (1) frc_data:   struct array with frc data to plot in the same graphs
%                 e.g. if you want to plot in the same graphs Sol1 and 
%                 Sol2, frc_data = {Sol1,Sol2}. (Sol1 and Sol2 must be in 
%                 the same format of output of 'nlvib_decode(..).m'. Extra
%                 field 'RMS' has to be included if requested in the plots.
% (2) items2plot: struct array with info on what to plot, in which figure 
%                 and in which format. It must have the following fields:
%                 (2.1) items2plot.subplot is a cell array indicating
%                 the suplots in figure. 
%                 e.g. items2plot.subplot = {221, 222, 223, 224}
%                 (2.2) items2plot.object is a cell array containing what
%                 to plot in each subplot. Cell arrays contains char
%                 vector: possible options are:
%                 'mod','phase','A','B','RMS','A0'
%                 e.g. items2plot.object = {'mod', 'mod', 'phase', 'phase'}
%                 (2.3) items2plot.dof is a vector containing the number of
%                 dofs that you want to plot in items2plot.subplot
%                 e.g. items2plot.dof = [1 2 1 2]
%                 (2.4) items2plot.harm is a vector array specifying which
%                 harmonics to plot is the subplots.
% (3) plotOptions:strucure array containing the following fields:
%                 (3.1)plotOptions.legend = cell array with legend for each
%                 subplot (same size of frc_data)
%                 (3.2) plotOptions.color = cell array with color for each
%                 curve plotted (same size of frc_data)
%                 (3.3) plotOptions.lineWidth = cell array with linewidth
%                 for each curve plotted (same size of frc_data)
%                 (3.4) plotOptions.lineStyle = cell array with linestyle
%                 for each curve plotted (same size of frc_data)
%                 (3.5) plotOptions.fontSize = fontsize for axes 
%                 (3.6) plotOptions.marker = cell array with specified
%                 marker for each curve plotted (same size of frc_curve).
%
% Author: Alexander Saccani, Msc in mechanical engineering
% University: Politecnico di Milano
% Created:  09/2021

%rename plot options
PlotColor = plotOptions.color;
PlotLineWidth = plotOptions.lineWidth;
PlotLineStyle = plotOptions.lineStyle;
PlotFontSize = plotOptions.fontSize;
PlotMarker = plotOptions.marker;
titleflag = isfield(plotOptions,'title');

%create figure
%figure('windowstyle','docked')

for n_subplot = 1:length(items2plot.subplot)
    
    %create subplot
    subplot(items2plot.subplot(n_subplot))

    for n_graph = 1:length(frc_data) %iterate on all data you want to plot

        frc = frc_data{n_graph};

        switch items2plot.object{n_subplot}

            case('A')

                dof2plot = items2plot.dof(n_subplot);
                harm = items2plot.harm(n_subplot);
                plot(frc.omega,frc.Qre(dof2plot,:,harm),'LineWidth',PlotLineWidth{n_graph},...
                        'LineStyle',PlotLineStyle{n_graph},'Marker',...
                        PlotMarker{n_graph},'color',PlotColor{n_graph});
                xlabel('\omega [rad/s]');
                ylabel('A');
                title(['A, dof ',num2str(dof2plot),', h=',num2str(harm)]);
                hold on

            case('B')

                dof2plot = items2plot.dof(n_subplot);
                harm = items2plot.harm(n_subplot);
                plot(frc.omega,frc.Qim(dof2plot,:,harm),'LineWidth',PlotLineWidth{n_graph},...
                        'LineStyle',PlotLineStyle{n_graph},'Marker',...
                        PlotMarker{n_graph},'color',PlotColor{n_graph});
                xlabel('\omega [rad/s]');
                ylabel('B ');
                title(['B, dof ',num2str(dof2plot),', h=',num2str(harm)]);
                hold on


            case('mod')

                dof2plot = items2plot.dof(n_subplot);
                harm = items2plot.harm(n_subplot);
                mod = abs(frc.Qre(dof2plot,:,harm) + 1i*frc.Qim(dof2plot,:,harm));
                plot(frc.omega,mod,'LineWidth',PlotLineWidth{n_graph},...
                        'LineStyle',PlotLineStyle{n_graph},'Marker',...
                        PlotMarker{n_graph},'color',PlotColor{n_graph});
                xlabel('\omega [rad/s]');
                ylabel('|A + iB|');
                title(['dof ',num2str(dof2plot),', h=',num2str(harm)]);
                hold on

            case('phase')

                dof2plot = items2plot.dof(n_subplot);
                harm = items2plot.harm(n_subplot);
                phase = angle(frc.Qre(dof2plot,:,harm) + 1i*frc.Qim(dof2plot,:,harm));
                plot(frc.omega,phase,'LineWidth',PlotLineWidth{n_graph},...
                        'LineStyle',PlotLineStyle{n_graph},'Marker',...
                        PlotMarker{n_graph},'color',PlotColor{n_graph});
                xlabel('\omega [rad/s]');
                ylabel('\angle(A + iB) [rad]');
                title(['dof ',num2str(dof2plot),', h=',num2str(harm)]);
                hold on

            case('RMS')

                dof2plot = items2plot.dof(n_subplot);
                plot(frc.omega,frc.RMS(dof2plot,:),'LineWidth',PlotLineWidth{n_graph},...
                        'LineStyle',PlotLineStyle{n_graph},'Marker',...
                        PlotMarker{n_graph},'color',PlotColor{n_graph});
                xlabel('\omega [rad/s]')
                ylabel('RMS')
                title(['dof ',num2str(dof2plot)]);
                hold on
                
            case('A0')
                
                dof2plot = items2plot.dof(n_subplot);
                plot(frc.omega,frc.Q0(dof2plot,:),'LineWidth',PlotLineWidth{n_graph},...
                        'LineStyle',PlotLineStyle{n_graph},'Marker',...
                        PlotMarker{n_graph},'color',PlotColor{n_graph});
                xlabel('\omega [rad/s]');
                ylabel('A_0');
                title(['A_0, dof ',num2str(dof2plot)]);
                hold on
                

            otherwise
                error('non valid object to plot, possible options are: A0,A,B,mod,phase,RMS')
        end

    end
    
    ax = gca;
    ax.FontSize = PlotFontSize;
    legend(plotOptions.legend, 'location','best');
    grid on
    if titleflag == 1
        title(plotOptions.title{n_subplot})
    end
    
end

end

        






















