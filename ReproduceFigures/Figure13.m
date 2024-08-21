clear all;
close all;

path(path,'../Models/Spiral_Wave/');
path(path,'Output_SpiralWave');

agrid = linspace(2.3,3.9,6);
bgrid = linspace(0.17,0.27,6);

U = cell(6,6);
% bgrid = bgrid(end:-1:1);
L = 250;
da = agrid(2) - agrid(1);
db = bgrid(2) - bgrid(1);

figure;
for i = 1:length(agrid)
    for j = 1:length(bgrid)
        
        cnew = num2str(agrid(i),'%1.6f');
        anew = num2str(bgrid(j),'%1.6f');
        fid = fopen('../Models/Spiral_Wave/Rossler/task.txt');
        C=textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        C{1}{3} = cnew; C{1}{1} = anew;
        writecell(C{1},'../Models/Spiral_Wave/Rossler/task.txt','QuoteStrings',0)
        
        cd ../Models/Spiral_Wave/Rossler
       
        system('./ezspiral')

        
        sol = readtable('fc.txt');
        sol = sol.Var2; 
        cd ../../../ReproduceFigures
        U{i,j} = reshape(sol(5:end),526,526);
        X = linspace(agrid(i) - 0.5*da, agrid(i) + 0.5*da, 526) ;
        Y = linspace(bgrid(j) - 0.5*db, bgrid(j) + 0.5*db, 526);
        imagesc(X,Y,U{i,j},'AlphaData', .5)
        hold on;
        
        xlim([agrid(1)-0.5*da,agrid(end)+0.5*da]);
        ylim([bgrid(1)-0.5*db,bgrid(end)+0.5*db]);
%         set(gcf,'position',[-130 1021 1176 796]);
        
        set(gca,'fontsize',24)
        
        
        xx = caxis; 
        if xx(2)-xx(1) < 1E-2
            caxis([xx(1)-0.1,xx(1)+0.1]);
        end
        set(gca,'YDir','normal')
        xlabel('$c$','interpreter','latex')
        ylabel('$a$','interpreter','latex')
        axis square
        pause(0.1)
        
        
    end
end
set(gcf,'position', [802 276 512 448])

load('result_Rossler_area_interp_bisect.mat')
plot(p_history{1}(1,:),p_history{1}(2,:),'^-','linewidth',2)
