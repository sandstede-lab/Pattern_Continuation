% main function for bifurcation tracing in Barkley model

clear all;
close all;
rng(5);

path(path,'../DataGenerator/Spiral_Wave/');

adapt = 1;

model = 'Barkley';
save_flag = 0;
steady = 1;
max_bisect = 2;

interp = 1;
dist = 'meander';
dist = 'retract';
sets = 'pos';
            
switch dist
    case 'meander'
    % line 1
        L = 0.5*2;


        init1 = [0.66,0.01]%+0.02;
        init2 = [0.67,0.02]%+0.02;

        max_step_arc_length = 0.01;
        min_step_arc_length = 0.01;
    % line 2
        L = 0.5*2;
        init1 = [0.39,0.01];
        init2 = [0.40,0.02];

        max_step_arc_length = 0.01;
        min_step_arc_length = 0.01;
        
    case 'drift'
        L = 0.5;
        init1 = [0.45,0.005];
        init2 = [0.48,0.02];
        max_step_arc_length = 0.005*1.5;
        min_step_arc_length = 0.005*1.5;

    case 'retract'
        L = 0.5;
        init1 = [0.37,0.01];
        init2 = [0.38,0.02];

        max_step_arc_length = 0.01;
        min_step_arc_length = 0.01;
end


[p_history,counts,L_history,p_history_all,metric_history_all] = bifurcation_tracing_bisection(init1, init2, L,max_step_arc_length, min_step_arc_length, dist, max_bisect, model, sets,steady,interp);


figure(1);
plot(p_history(1,:), p_history(2,:), 'o-' , 'linewidth',2);
legend()
hold on;

title(dist)

ylim([0.01,0.13])
xlim([0.25,0.9])
xlabel('a')
ylabel('b')


set(gca,'fontsize',24)




if save_flag

    save(['result_',model,'_',dist,'_adaptsize_interp_bisect.mat'],'p_history','counts','dist','p_history_all','L_history','metric_history_all');




end




