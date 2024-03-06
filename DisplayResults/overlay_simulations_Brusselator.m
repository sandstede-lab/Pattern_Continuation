clear all;
close all;

% Reproduce the figures in paper (result+simulated patterns)

path(path,'../DataGenerator/Reaction_Diffusion');

agrid = linspace(3.5,7,7);
bgrid = linspace(6,14,7);
bgrid = bgrid(end:-1:1);
da = agrid(2) - agrid(1);
db = bgrid(1) - bgrid(2);


show_sim = 1;
    
figure;
for i = 1:length(agrid)
    for j = 1:length(bgrid)
        if show_sim
%             U{i,j} = reaction_diffusion_eqn_brusselator_2(agrid(i), bgrid(j));
            U{i,j} = Brusselator(agrid(i), bgrid(j), 0, 50, 200, agrid(i)/2+5+1*randn(50), bgrid(i)/agrid(i)/2+5+1*randn(50));
            
            
            if max( U{i,j}(:) ) - min( U{i,j}(:)) > 0.01
                U{i,j} = reshape( zscore(U{i,j}(:)), size(U{i,j},1),size(U{i,j},2) );
            else
                U{i,j} = 0*U{i,j};
            end

            X = linspace(agrid(i) - 0.5*da, agrid(i) + 0.5*da, size(U{i,j},1)) ;
            Y = linspace(bgrid(j) - 0.5*db, bgrid(j) + 0.5*db, size(U{i,j},2)) ;
            
%             X = linspace(agrid(i) - 0.48*da, agrid(i) + 0.48*da, size(U{i,j},1)) ;
%             Y = linspace(bgrid(j) - 0.48*db, bgrid(j) + 0.48*db, size(U{i,j},2)) ;


            imagesc(X,Y,U{i,j}, 'AlphaData', .5); hold on;
        end
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(end)-0.5*db,bgrid(1)+0.5*db]);
        
        set(gca,'fontsize',16)
        
        
        xx = caxis; 
        if xx(2)-xx(1) < 1E-2
            caxis([xx(1)-0.1,xx(1)+0.1]);
        end
        set(gca,'YDir','normal')
        
        axis square
        pause(0.1)
        xlabel('$a$','interpreter','latex')
        ylabel('$b$','interpreter','latex')
    end
end

C= brewermap([],'RdYlBu'); colormap(C);
caxis([-1.5,1.5])


set(gcf,'position', [802 276 512 448])
load('result_Brusselator_interp_pos_steady_bisection2.mat')
plot(p_history{1}(1,:),p_history{1}(2,:),'ko-','linewidth',2)
load('result_Brusselator_interp_pos_bisection2.mat')
plot(p_history{2}(1,:),p_history{2}(2,:),'o-','linewidth',2)
load('result_Brusselator_interp_neg_bisection2.mat')
plot(p_history{2}(1,:),p_history{2}(2,:),'^-','linewidth',2)


if show_sim == 0
    figure;
    for i = 1:2
        figure(i)
        load('result_Brusselator_interp_pos_steady_bisection2.mat')
        plot(p_history{1}(1,:),p_history{1}(2,:),'ko-','linewidth',2),hold on,
        load('result_Brusselator_interp_pos_bisection2.mat')
        plot(p_history{i}(1,:),p_history{i}(2,:),'bo-','linewidth',2)
        load('result_Brusselator_interp_neg_bisection2.mat')
        plot(p_history{i}(1,:),p_history{i}(2,:),'^-','linewidth',2)
        load('result_Brusselator_interp_pos_steady_bisection2.mat')
        scatter(p_history_all{1}(1,:),p_history_all{1}(2,:),25,metric_history_all{1},'*')

        load('result_Brusselator_interp_pos_bisection2.mat')
        scatter(p_history_all{i}(1,:),p_history_all{i}(2,:),25,metric_history_all{i},'o','filled')
        load('result_Brusselator_interp_neg_bisection2.mat')
        scatter(p_history_all{i}(1,:),p_history_all{i}(2,:),25,metric_history_all{i},'^','filled')
        C= brewermap([],'RdYlBu'); colormap(C);
%         caxis([-1.5,1.5])
        colorbar()
        if max( metric_history_all{i} ) > 1
            title('Metric: number of $\alpha$-shapes','interpreter','latex')
        else
            title('Metric: roundness of $\alpha$-shapes','interpreter','latex')
        end
        set(gcf,'position', [802 276 512 448])
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(end)-0.5*db,bgrid(1)+0.5*db]);
        set(gca,'fontsize',24)
        grid on
        axis square
        pause(0.1)
        xlabel('$a$','interpreter','latex')
        ylabel('$b$','interpreter','latex')
    end
    
    
    
end



   