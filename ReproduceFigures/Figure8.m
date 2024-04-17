clear all;
close all;

path(path,'../DataGenerator/Reaction_Diffusion');
path(path,'Output_ReactionDiffusion2D');




pattern_init = 'spots';

        
agrid = linspace(-0.6,0.2,7);
bgrid = linspace(0.2,3,7);

bgrid = bgrid(end:-1:1);
da = agrid(2) - agrid(1);
db = bgrid(1) - bgrid(2);
        
        
figure;
for i = 1:length(agrid)
    for j = 1:length(bgrid)
        
        U{i,j} = SH_2D(agrid(i), bgrid(j),pattern_init);
        N = size(U{i,j},1);
        U{i,j} = U{i,j}(1:N/2,1:N/2);
        if max(U{i,j}(:)) - min(U{i,j}(:)) > 0.01
            U{i,j} = reshape(zscore(U{i,j}(:)),N/2,N/2);
        else
            U{i,j} = 0*U{i,j};
        end
        
        X = linspace(agrid(i) - 0.5*da, agrid(i) + 0.5*da, size(U{i,j},1)) ;
        Y = linspace(bgrid(j) - 0.5*db, bgrid(j) + 0.5*db, size(U{i,j},2)) ;
        
        
        
        imagesc(X,Y,U{i,j}, 'AlphaData', .5); hold on;
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(end)-0.5*db,bgrid(1)+0.5*db]);
%         set(gcf,'position',[-130 1021 1176 796]);
        
        set(gca,'fontsize',24)
        
        
        xx = caxis; 
        if xx(2)-xx(1) < 1E-2
            caxis([xx(1)-0.1,xx(1)+0.1]);
        end
        set(gca,'YDir','normal')
        xlabel('$\mu$','interpreter','latex')
        ylabel('$\nu$','interpreter','latex')
        axis square
        pause(0.1)
        
        
    end
end

C = brewermap([],'RdYlBu'); colormap(C);
caxis([-1.5,1.5]);
set(gcf,'position', [802 276 512 448])
load(['result_SH_',pattern_init,'_all_bisection.mat'])
plot(p_history{1}(1,:),p_history{1}(2,:),'ko-','linewidth',2)





