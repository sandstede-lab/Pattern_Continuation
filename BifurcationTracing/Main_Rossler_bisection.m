% main function for bifurcation tracing in Rossler

clear all;
close all;
rng(5);


path(path,'../DataGenerator/Spiral_Wave/');

adapt = 1;

model = 'Rossler';
save_flag = 0;
steady = 0;
max_bisect = 2;

interp = 1;
dist = 'area';
sets = 'pos';
                
                
L = 1.5;
init1 = [3.9,0.17];
init2 = [3.8,0.17];

max_step_arc_length = 0.03;
min_step_arc_length = 0.03;

[p_history,counts,L_history,p_history_all,metric_history_all] = bifurcation_tracing_bisection(init1, init2, L,max_step_arc_length, min_step_arc_length, dist, max_bisect, model, sets,steady, interp);

                   

figure(1);
plot(p_history(1,:), p_history(2,:), linestyles , 'linewidth',2);
legend()
hold on;

title(dist)

xlim([2.5,4])
ylim([0.16,0.26])
xlabel('c')
ylabel('a')


set(gca,'fontsize',24)





if save_flag
    save(['result_',model,'_',dist,'_interp_bisect.mat'],'p_history','counts','dists','p_history_all','L_history','metric_history_all');


end





  
        
  
          