clear all;
close all;

path(path, '../Models/Reaction_Diffusion')
path(path,'Output_ReactionDiffusion2D');

figure;

load('result_Brusselator_interp_pos_bisection2.mat')

subplot(2,4,[2,3,6,7]);
plot(p_history{1}(1,:),p_history{1}(2,:),'k--','linewidth',2)
xlim([3,6])
ylim([6,15])
xticks([]);
yticks([]);
xlabel('$a$','interpreter','latex');
ylabel('$b$','interpreter','latex')

hold on,

scatter(p_history_all{1}(1,:),p_history_all{1}(2,:),[],metric_history_all{1},'filled');

scatter(4.02195,10.5687,'*')
scatter(5.06292,10.7525,'*')
set(gca,'fontsize',24)
caxis([0,1])
colorbar()

subplot(2,4,8)
% a = 4.5070; b = 8.6742;
a = 5.06292; b = 10.7525;

U = Brusselator(a,b);
imagesc(U);
xticks([]);
yticks([]);
caxis([-3,11])
cm = brewermap([],'RdYlBu');
colormap(cm)

roundness = [];
Nx = 50;
[X,Y] = meshgrid( linspace(0,50,Nx), linspace(0,50,Nx));
X = X(:); Y = Y(:);
alpha = 1;

for j = 1:10
    Ubp{j} = Brusselator(a,b);
    idx = find(Ubp{j} > quantile(Ubp{j}(:), 0.7));
    shp_pos = alphaShape(X(idx),Y(idx),alpha);
    areas = area(shp_pos, 1:numRegions(shp_pos));
    perimeters = perimeter(shp_pos, 1:numRegions(shp_pos));   
    roundness = [roundness, 4*pi*areas./(perimeters.^2)];
    
    
end



subplot(2,4,4)
histogram(roundness,'BinEdges',linspace(0,1,20))
title({'histogram: roundness',[num2str(length(roundness)),' $\alpha$-shapes']},'interpreter','latex');
set(gca,'fontsize',16)
xlim([0,1])

subplot(2,4,1)
% a = 3.9257; b = 9.6796;
a = 4.02195; b = 10.5687;
U = Brusselator(a,b);
imagesc(U);
xticks([]);
yticks([]);
caxis([-2,9])
colormap(cm)

roundness = [];
Nx = 50;
[X,Y] = meshgrid( linspace(0,50,Nx), linspace(0,50,Nx));
X = X(:); Y = Y(:);
alpha = 1;

for j = 1:10
    Ubp{j} = Brusselator(a,b);
    idx = find(Ubp{j} > quantile(Ubp{j}(:), 0.7));
    shp_pos = alphaShape(X(idx),Y(idx),alpha);
    areas = area(shp_pos, 1:numRegions(shp_pos));
    perimeters = perimeter(shp_pos, 1:numRegions(shp_pos));   
    roundness = [roundness, 4*pi*areas./(perimeters.^2)];
    
    
end

subplot(2,4,5)
histogram(roundness,'BinEdges',linspace(0,1,20))
title({'histogram: roundness',[num2str(length(roundness)),' $\alpha$-shapes']},'interpreter','latex');
set(gca,'fontsize',16)
xlim([0,1])

set(gcf,'position',[24 356 889 371])
