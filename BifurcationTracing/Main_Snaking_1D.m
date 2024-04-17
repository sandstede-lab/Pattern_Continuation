% This file shows code for snaking in 1D Swift-Hohenberg

clear all;
close all;
rng(5);

path(path,'../DataGenerator/Reaction_Diffusion');

dist = 'num';% distance metric to be used. 

interp = 1; % bifurcation tracing with interpolation to improve smoothness
adapt = 1; % arc length per step is adaptive

sets = 'pos';

save_flag = 0; % results will be saved

model = 'SH_1D';

steady = 0; % trace out when homogeneous solution occurs

max_bisect = 2; % Number of bisection allowed
        
    

max_step_arc_length = 0.02;
min_step_arc_length = max_step_arc_length;
L = 1; % Arc length in total

if steady == 0
    init1 = [0.3,1.9]; init2 = [0.3,2];
else

    init1 = [0.36,1.7]; init2 = [0.36,1.8];
end
    
    
  
                    
[p_history,counts,L_history,p_history_all,metric_history_all] = bifurcation_tracing_bisection(init1, init2, L, max_step_arc_length, min_step_arc_length, dist,max_bisect,model,sets,steady,interp);



figure;
%                 figure(j);
plot(p_history(1,:), p_history(2,:), 'o-' , 'linewidth',2);
legend()
hold on;

                


ylim([1.8,2.7])
xlim([0,1.2])
xlabel('$\mu$','interpreter','latex')
ylabel('$\nu$','interpreter','latex')

    
set(gca,'fontsize',24)
               
                
%% Save results as a file if necessary

        
if save_flag
   if steady

       save(['result_',model,'_interp_',sets,'_steady_bisection',num2str(max_bisect),'.mat'],'p_history','counts','dist','metric_history_all','L_history','p_history_all');
   else
       save(['result_',model,'_interp_',sets,'_bisection',num2str(max_bisect),'.mat'],'p_history','counts','dist','metric_history_all','L_history','p_history_all');

   end


end
    