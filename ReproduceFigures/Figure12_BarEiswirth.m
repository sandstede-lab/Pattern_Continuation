clear all;
close all;

path(path,'../DataGenerator/Spiral_Wave/');
path(path,'Output_SpiralWave');

bgrid = linspace(0.04,0.24,2);
agrid = linspace(0.01,0.08,2);
tipping_pts = cell(length(agrid),length(bgrid));
% bgrid = bgrid(end:-1:1);
L = 50;
da = agrid(2) - agrid(1);
db = bgrid(2) - bgrid(1);

figure;
for i = 1:length(agrid)
    for j = 1:length(bgrid)
        
        epnew = num2str(1/agrid(i),'%.6f');
        bnew = num2str(bgrid(j),'%1.6f');
        fid = fopen('../DataGenerator/Spiral_Wave/Bar-Eiswirth/task.txt');
        C=textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        C{1}{3} = epnew; C{1}{2} = bnew;
        writecell(C{1},'../DataGenerator/Spiral_Wave/Bar-Eiswirth/task.txt','QuoteStrings',0)
        
        cd ../DataGenerator/Spiral_Wave/Bar-Eiswirth
        system('rm tip.txt')
        system('touch tip.txt')
        system('./ezspiral')


        cd ../../../ReproduceFigures
        tipping_pts{i,j} = readtable('../DataGenerator/Spiral_Wave/Bar-Eiswirth/tip.txt');
        if size(tipping_pts{i,j})
            tipping_pts{i,j} = [tipping_pts{i,j}.Var2';tipping_pts{i,j}.Var3'];
            tipping_pts{i,j} = tipping_pts{i,j}(:,501:end);
        end
        
        X = [agrid(i) - 0.5*da, agrid(i) + 0.5*da];
        Y = [bgrid(j) - 0.5*db, bgrid(j) + 0.5*db];
        
        if prod(size(tipping_pts{i,j})) > 0
            
            tipping_pts{i,j}(1,:) = (0.5*(tipping_pts{i,j}(1,:)-min(tipping_pts{i,j}(1,:)))/( max(tipping_pts{i,j}(1,:))-min(tipping_pts{i,j}(1,:)))+0.25)*(X(end)-X(1))+X(1);
            tipping_pts{i,j}(2,:) = (0.5*(tipping_pts{i,j}(2,:)-min(tipping_pts{i,j}(2,:)))/( max(tipping_pts{i,j}(2,:))-min(tipping_pts{i,j}(2,:)))+0.25)*(Y(end)-Y(1))+Y(1);
            
            scatter(tipping_pts{i,j}(1,:),tipping_pts{i,j}(2,:),6,'b','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);hold on;
        end
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(1)-0.5*db,bgrid(end)+0.5*db]);
%         set(gcf,'position',[-130 1021 1176 796]);
        
        set(gca,'fontsize',24)
        
        
        xx = caxis; 
        if xx(2)-xx(1) < 1E-2
            caxis([xx(1)-0.1,xx(1)+0.1]);
        end
        set(gca,'YDir','normal')
        xlabel('$\epsilon$','interpreter','latex')
        ylabel('$b$','interpreter','latex')
        axis square
        pause(0.1)
        
        
    end
end
set(gcf,'position', [802 276 512 448])

% set(gcf,'position', [802 276 512 448])
load('result_BE_ring_fit_adaptsize_interp_bisect.mat')
for i = 2:3
    plot(p_history{i}(1,:),p_history{i}(2,:),'o-','linewidth',2)
end

load('result_BE_turbulence2_adaptsize_interp_bisect.mat')
plot(p_history{1}(1,:),p_history{1}(2,:),'o-','linewidth',2)
load('result_BE_adaptsize_steady_interp_bisect.mat')
plot(p_history{1}(1,:),p_history{1}(2,:),'ko-','linewidth',2)