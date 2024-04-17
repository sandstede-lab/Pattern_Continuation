clear all;
close all;

warning off

path(path,'Output_Freezing');


radial = 1;
domain_size = 30;

params = 'm';

switch params
    case 'd'
        load results_d
    case 'm'
        load results_m
    case 'r'
        load results_r
end

% total_steps = 1200000;
% 
step = 4000; % align spirals every 4000 steps
% r = 20; % radius of area to be matched
    

figure;
subplot(1,2,1)
theta_c = cumsum(real(theta));
plot(0.05*step*[1:length(theta_c)],theta_c,'linewidth',2)

title('$\theta$','interpreter','latex')
grid on
xlabel('$t$','interpreter','latex')
set(gca,'fontsize',24)

subplot(1,2,2)
for t = 1:size(b,2)
    
    R = [cos(theta(t)), sin(theta(t)); -sin(theta(t)), cos(theta(t))];

    bt(:,t) = R*b(:,t);
    for tt = t-1:-1:1
        R = [cos(theta(tt)), sin(theta(tt)); -sin(theta(tt)), cos(theta(tt))];
        bt(:,t) = R* (bt(:,t) + b(:,tt));
    end
    

    
end

dx = diff(bt(1,:)); dy = diff(bt(2,:));
quiver(bt(1,1:end-1),bt(2,1:end-1),dx,dy,0,':','linewidth',2,'color',[0.7,0.7,0.7]), hold on,
scatter(bt(1,:),bt(2,:), [], 1:size(bt,2),'filled');
title('$\vec b(t)$','interpreter','latex')
xlabel('$b_1(t)$','interpreter','latex')
ylabel('$b_2(t)$','interpreter','latex')
set(gca,'fontsize',24)
grid on
set(gcf,'position',[440 434 874 363])







