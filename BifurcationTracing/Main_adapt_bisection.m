% main function for bifurcation tracing

clear all;
close all;
rng(5);


dist = 'num';% distance metric to be used. Alternative: 'roundness'

interp = 1; % bifurcation tracing with interpolation to improve smoothness
adapt = 1; % arc length per step is adaptive

sets = 'pos';

save_flag = 0; % results will be saved

model = 'Brusselator'; % name of IVP solver. Alternative: 'SH','Gray_Scott','Schnakenberg'

steady = 1; % trace out when homogeneous solution occurs

max_bisect = 2; % Number of bisection allowed
        
    
switch model
    case 'Brusselator'
        max_step_arc_length = 0.5;
        min_step_arc_length = max_step_arc_length/2;

        switch sets
            case 'pos'
                init1 = [3.5,6]; init2 = [3.5,6.5];

            case 'neg'

                init1 = [4.2,6]; init2 = [4.2,7];
        end
        if steady
            % init1 = [4.5,6]; init2 = [4.5,6.5]; % Trace from
            % below
            init1 = [7.8,14.5]; init2 = [7.8,14];
        end
        L = 8;
    case 'SH'
        max_step_arc_length = 0.15;
        min_step_arc_length = max_step_arc_length/2;
        L = 1.8*1.5;
        if steady == 0
            init1 = [0.0,0.15]; init2 = [0.02,0.15];
        else
            max_step_arc_length = 0.15;
            min_step_arc_length = max_step_arc_length/2;
            init1 = [0.0,0]; init2 = [0.0,0.1];
            init1 = [0.0,0.3]; init2 = [0.0,0.4];
        end
    

    case 'GS'
        max_step_arc_length = 0.0012;
        min_step_arc_length = max_step_arc_length/2;
        switch sets
            case 'neg'
                init1 = [0.056,0.024]; init2 = [0.056,0.025];
            case 'pos'
                init1 = [0.053,0.024]; init2 = [0.053,0.025];
        end
        L = 0.018;
    
    case 'Schnakenberg'
        max_step_arc_length = 0.5;
        min_step_arc_length = max_step_arc_length/2;
        switch sets
            case 'neg'
                init1 =[6.5,6]; init2 = [6.6,6];
            case 'pos'
                
                init1 =[4.9,4.8]; init2 = [5,4.8];
        end
        if steady
             init1 =[4.3,5.3]; init2 = [4.5,5.3];
        end

        L = 6;

end

linestyles = {'o-', '*--', 'x:','^-.'};

  
            
  
                    
[p_history,counts,L_history,p_history_all,metric_history_all] = bifurcation_tracing_bisection(init1, init2, L, max_step_arc_length, min_step_arc_length, dist,max_bisect,model,sets,steady,interp);



figure;
%                 figure(j);
plot(p_history(1,:), p_history(2,:), linestyles , 'linewidth',2);
legend()
hold on;

                

switch model
    case 'Brusselator'

        ylim([6,13])
        xlim([3,7])
        xlabel('a')
        ylabel('b')
    case 'SH'
        ylim([0,2])
        xlim([-0.15,0.4])
        xlabel('\mu')
        ylabel('\nu')

    case 'GS'
        xlim([0.05,0.065])
        ylim([0.023,0.038])
        ylabel('F')
        xlabel('k')


    case 'Schnakenberg'
        xlim([4.5,7.5]);
        ylim([4.5,6.5]);
        xlabel('a')
        ylabel('b')
end
set(gca,'fontsize',24)
               
                
%% Save results as a file if necessary

        
if save_flag
   if steady

       save(['result_',model,'_interp_',sets,'_steady_bisection',num2str(max_bisect),'.mat'],'p_history','counts','dist','metric_history_all','L_history','p_history_all');
   else
       save(['result_',model,'_interp_',sets,'_bisection',num2str(max_bisect),'.mat'],'p_history','counts','dist','metric_history_all','L_history','p_history_all');

   end


end
    