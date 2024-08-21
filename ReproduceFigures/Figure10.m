clear all;
close all;

path(path, '../Models/Reaction_Diffusion')
path(path,'Output_Snaking');

agrid = linspace(0.3,1.1,7);
bgrid = linspace(1.7,2.6,7);

bgrid = bgrid(end:-1:1);
da = agrid(2) - agrid(1);
db = bgrid(1) - bgrid(2);


show_sim = 1;
    
figure;
for i = 1:length(agrid)
    for j = 1:length(bgrid)
        if show_sim
            U{i,j} = SH_1D(agrid(i), bgrid(j));
            

            X = linspace(agrid(i) - 0.5*da, agrid(i) + 0.5*da, size(U{i,j},1)) ;
            Y = linspace(bgrid(j) - 0.5*db, bgrid(j) + 0.5*db, size(U{i,j},1)) ;
            if max(U{i,j}) - min(U{i,j}) > 0.001
                U{i,j} = (0.5*(U{i,j}-min(U{i,j}))/( max(U{i,j})-min(U{i,j}))+0.25)*(Y(end)-Y(1))+Y(1);
            else
                U{i,j} = 0*U{i,j};
            end
            U{i,j} = U{i,j} - U{i,j}(1) + 0.5*(Y(end)-Y(1))+Y(1);
            
            plot(X,U{i,j},'Color',[0.2 0.5 0.9 0.5]);hold on;
        
            
        end
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(end)-0.5*db,bgrid(1)+0.5*db]);
%         set(gcf,'position',[-130 1021 1176 796]);
        
        set(gca,'fontsize',16)
        
        
        xx = caxis; 
        if xx(2)-xx(1) < 1E-2
            caxis([xx(1)-0.1,xx(1)+0.1]);
        end
        set(gca,'YDir','normal')
        
        axis square
        pause(0.1)
        xlabel('$\mu$','interpreter','latex')
        ylabel('$\nu$','interpreter','latex')
    end
end

C= brewermap([],'RdYlBu'); colormap(C);
caxis([-1.5,1.5])


set(gcf,'position', [802 276 512 448])
load('result_SH_adaptsize_pos_steady_bisect.mat')
plot(p_history{1}(1,:),p_history{1}(2,:),'k.-','markersize',12,'linewidth',2),hold on
load('result_SH_adaptsize_pos_bisect.mat')
plot(p_history{1}(1,:),p_history{1}(2,:),'.-','markersize',12,'linewidth',2)

AUTO1 = readtable( 'Snaking_Diagram_1n.txt');
hold on,
plot(AUTO1.Var5, AUTO1.Var11,'color',[.5,.5,.5],'linewidth',2);
AUTO2 = readtable('Snaking_Diagram_2.txt');
hold on,
plot(AUTO2.Var5, AUTO2.Var11,'color',[.5,.5,.5],'linewidth',2);



% figure;
if show_sim == 0
    for i = 1
        figure(i)
        load('result_SH_adaptsize_pos_steady_bisect.mat')
        plot(p_history{1}(1,:),p_history{1}(2,:),'ko-','linewidth',2),hold on,
        load('result_SH_adaptsize_pos_bisect.mat')
        plot(p_history{i}(1,:),p_history{i}(2,:),'bo-','linewidth',2)
       

        load('result_SH_adaptsize_pos_steady_bisect.mat')
        scatter(p_history_all{i}(1,:),p_history_all{i}(2,:),25,metric_history_all{i},'o','filled')
        load('result_SH_adaptsize_pos_bisect.mat')
        scatter(p_history_all{i}(1,:),p_history_all{i}(2,:),25,metric_history_all{i},'^','filled')
        C= brewermap([],'RdYlBu'); colormap(C);
%         caxis([-1.5,1.5])
        colorbar()
        
        set(gcf,'position', [802 276 512 448])
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(end)-0.5*db,bgrid(1)+0.5*db]);
        set(gca,'fontsize',24)
        grid on
        axis square
        pause(0.1)
        xlabel('$\mu$','interpreter','latex')
        ylabel('$\nu$','interpreter','latex')
    end
    
    
    
end



   