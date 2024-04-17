% main function for bifurcation tracing

clear all;
close all;
rng(5);


interp = 1;
adapt = 1;
% models = {'SH','Brusselator'};
sets = 'pos';
model = 'Bar-Eiswirth';%,'SH'}

i = 1;
save_flag = 0

steady = 1;

max_bisect = 2;
interp = 1;




switch steady
    case 0

       dists = {'meander','meander2'}%,'meander'}
%                dists = {'meander2'};
       dists = {'turbulence','meander'};
       dists = {'turbulence'};
    case 1
        dists = {'retract'}
end

    linestyles = {'o-', '*--', 'x:','^-.'};

    for j = 1:length(dists)

            switch dists{j}
                case 'meander'
                    L = 0.5;

                        % line 1:
                        init1 = [0.058,0.04];
                        init2 = [0.058,0.05];


                        max_step_arc_length = 0.005;
                        min_step_arc_length = 0.005;
                        
                        % line 2:
                        
                        init1 = [0.006,0.1];
                        init2 = [0.006,0.11];


                        max_step_arc_length = 0.003;
                        min_step_arc_length = 0.003;
                        
                case 'retract'
                        L = 0.1;
                        init1 = [0.19,0.2];
                        init2 = [0.18,0.2];
%                             
                        init1 = [0.07,0.19];
                        init2 = [0.06,0.19];

                        max_step_arc_length = 0.005;
                        min_step_arc_length = 0.005;
                            
                        
                
       

                        
                case 'turbulence'

                    L = 0.5;


                    init1 = [0.07,0.03];
                    init2 = [0.07,0.04];




                    max_step_arc_length = 0.005;
                    min_step_arc_length = 0.005;


            end
                

            
            [p_history{j},counts{j},L_history{j},p_history_all{j},metric_history_all{j}] = bifurcation_tracing_bisection(init1, init2, L,max_step_arc_length, min_step_arc_length, dists{j}, max_bisect, model, sets,steady,interp);

                
            

%                 [p_history{j},counts{j},L_history{j},p_history_all{j},metric_history_all{j}] = bifurcation_tracing_bisection(init1, init2, L,max_step_arc_length, min_step_arc_length, dists{j}, max_bisect, model, sets,steady);


            figure(1);
%                 figure(j);
            plot(p_history{j}(1,:), p_history{j}(2,:), linestyles{j} , 'linewidth',2);
            legend()
            hold on;

            title(dists{j})

            switch model
                case 'Barkley'

                    ylim([0.01,0.13])
                    xlim([0.25,0.9])
                    xlabel('a')
                    ylabel('b')
                case 'Bar-Eiswirth'
                    xlim([0.0,0.15])
                    ylim([0.04,0.26])
                    xlabel('a')
                    ylabel('\epsilon')
                case 'Rossler'
                    xlim([2.5,4])
                    ylim([0.16,0.26])
                    xlabel('c')
                    ylabel('a')

            end
            set(gca,'fontsize',24)

            label = cell(1,j);
            for ii = 1:j
                label{ii} = dists{ii};
            end
            legend(label,'Location','northwest')


            pause(0.1);


           
        if save_flag
            if interp == 0
               if steady

                   save(['result_',model,'_adaptsize_steady_bisect.mat'],'p_history','counts','dists','p_history_all','L_history','metric_history_all');

               else
                   save(['result_',model,'_',dists{j},'_adaptsize_bisect.mat'],'p_history','counts','dists','p_history_all','L_history','metric_history_all');
               end
               
            else
                if steady

                   save(['result_',model,'_adaptsize_steady_interp_bisect.mat'],'p_history','counts','dists','p_history_all','L_history','metric_history_all');

               else
                   save(['result_',model,'_',dists{j},'_adaptsize_interp_bisect.mat'],'p_history','counts','dists','p_history_all','L_history','metric_history_all');
                end
            end
           

            
        end





    end


% save(['result_',model,'_adaptsize_',sets,'_',solver,'_line.mat'],'p_history','counts','dists','val_history','p_history_all','L_history','p_history_all_sim', 'val_feat_all_sim');
           