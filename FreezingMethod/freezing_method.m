clear all;
close all;

warning off

path(path,'../DataGenerator/Spiral_Wave/');
path(path,'.');

new_file = 1;
radial = 1;
domain_size = 30;

params = 'm';
make_movie = 0;
if make_movie
    movie_name = ['results_',params,'.mp4'];

    writerObj = VideoWriter(movie_name,'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
end

total_steps = 1200000;

step = 4000; % align spirals every 4000 steps
r = 20; % radius of area to be matched
    
U_history = cell(1, round(total_steps/step) );
tip_pt = zeros(2,  round(total_steps/step) );
counter = 1;
cd ../DataGenerator/Spiral_Wave/Barkley_radial_c
system('rm ic.dat');

fid = fopen('task.dat');
C=textscan(fid,'%s','delimiter','\n');
fclose(fid);



switch params
    case 'm'
        switch domain_size
            case 30
                system('cp fc_m.dat ic.dat')
            case 50
                system('cp fc_m_big.dat ic.dat')
        end

        
        C{1}{1} = 0.71; 
        
        
        
    case 'r'
        
        switch domain_size
            case 30
                system('cp fc_r.dat ic.dat')
            case 50
                system('cp fc_r_big.dat ic.dat')
        end

        
        C{1}{1} = 0.8; 
        C{1}{1} = 0.82; 
        
        

    case 'd'
        switch domain_size
            case 30
                system('cp fc_d.dat ic.dat')
            case 50
                system('cp fc_d_big.dat ic.dat')
        end
        
        C{1}{1} = 0.63; 
        
       
        
end
C{1}{12} = step;
C{1}{4} = domain_size;
writecell(C{1},'task.dat','QuoteStrings',0)
h = 1/2;


a = readtable('ic.dat');
aa = a.Var1;   
a = a.Var2;


% domain discretization
C{1}{8} = 161;
C{1}{9} = 160;
XX = linspace(0,2*pi,C{1}{9}+1); XX = XX(1:end-1); 
YY = linspace(0,30,C{1}{8}+1); YY = YY(2:end); 
X = YY.*cos(XX)';
Y = YY.*sin(XX)';

U0 = aa(5:end); V0 = a(5:end);
U_prior = U0; V_prior = V0;
U0 = reshape(U0(1:end-1),   size(X,1), size(X,2));
V0 = reshape(V0(1:end-1),   size(X,1), size(X,2));
tip0 = [0,0]; tip_current = [0,0];

theta = zeros(1,round(total_steps/step));
b = zeros(2,round(total_steps/step));
for i = 1:round(total_steps/step)
    i
    
    
    [U_current,V_current] = ezspiral_my(U_prior,V_prior);

    
    cd ../../../FreezingMethod

    
    
    % Figure out rotation/translation
    if i > 1
        init = [thetat; b(:,i-1)];
    else
        init = [0,0,0]';
    end

    [R,b(:,i),x1,x2,thetat] = find_rotation_translation_alpha( U0, U_current, [X(:);0], [Y(:);0], tip0, tip_current, r, h,init);
    

    cd ../DataGenerator/Spiral_Wave/Barkley_radial_c
    
    
    
%     bt(:,i) = inv(eye(2)-R)*b(:,i);
    theta(i) = acos(R(1,1));
    theta(i) = real(theta(i));
    if sign(-sin(theta(i))) ~= sign( R(1,2))
        theta(i) = -theta(i);
        
    end
    if i > 1 
        
        theta(i) = theta(i-1) - pi + mod(theta(i)-theta(i-1)-pi,2*pi);
    end
    figure(1);
    
    scatter(x1(:,1),x1(:,2),'b'), hold on,
    
    scatter(x2(:,1),x2(:,2),'r'), 
    x1_map = (b(:,i)' + x1) *R';
    
%     tip_map = b(:,i)' + tip0'*R';
    scatter( x1_map(:,1), x1_map(:,2),'y')

    hold off
    set(gca,'fontsize',12)
    
    % Modify intial condition according to the rotation/translation
    
    S_current = [[X(:);0],[Y(:);0]]*R - b(:,i)';
    
%     S_current = [X(:),Y(:)];
    X_current = S_current(:,1);
    Y_current = S_current(:,2);
      
            
    
    
    figure(2)
    subplot(1,4,1)
%     mesh(X,Y,U0)
%     set(gca,'view',[0,90])
    scatter(X(:),Y(:),[],U0(:),'filled');
    title('Solution at 0')
    set(gca,'fontsize',12)
%     scatter(X(:),Y(:),[],U0(:))
    
    subplot(1,4,2)
%     mesh(X,Y,U_current)
%     set(gca,'view',[0,90])
    scatter([X(:);0],[Y(:);0],[],U_current(:),'filled');
    title('Solution at t')
    set(gca,'fontsize',12)
%     scatter(X(:),Y(:),[],U_current(:));
    
    subplot(1,4,3)
%     mesh(X_current,Y_current,U_current)
%     set(gca,'view',[0,90])
    scatter(X_current(:),Y_current(:),[],U_current(:),'filled');
    title('Solution at t rotation')
    xlim([-30,30])
    ylim([-30,30])
    
    set(gca,'fontsize',12)
%     scatter(X_current(:),Y_current(:),[],U_current(:));
    
    subplot(1,4,4)

    old_U = U_current;
    old_V = V_current;
    

    
    new_U = griddata(S_current(:,1),S_current(:,2), U_current(:), [X(:);0],[Y(:);0]);
    new_V = griddata(S_current(:,1),S_current(:,2), V_current(:), [X(:);0],[Y(:);0]);
    idx = find(isnan(new_U));
    
    
    
    new_U(idx) = old_U(idx);
    new_V(idx) = old_V(idx);
    
    U_prior = new_U;
    V_prior = new_V;
    
%     mesh(X,Y,reshape(new_U,size(X,1),size(X,2)));
%     
%     set(gca,'view',[0,90])
    scatter([X(:);0],[Y(:);0],[],new_U,'filled');
    set(gca,'fontsize',12)
    title('Solution at t rotation+interpolation')

    set(gcf,'position',[40 512 1387 263])
    pause(0.1)
    
   

    
    
    
    
    U_history{i} = sparse(U_prior);
    
    if make_movie
        figure(3);
        subplot(3,2,1)
        mesh(X,Y,U0)
        set(gca,'view',[0,90])
        title('Solution at 0')
        set(gca,'fontsize',12)
        
        subplot(3,2,2)
        mesh(X,Y,reshape(new_U(1:end-1), size(X,1), size(X,2)));
        set(gca,'view',[0,90])
        set(gca,'fontsize',12)
        title('Solution at t')
    
        subplot(3,2,3)
        %plot(theta(1:i),'linewidth',2)
        
        plot(b(1,1:i),'linewidth',2) , hold on, plot(b(2,1:i),'linewidth',2)
        
        hold off
        xlabel('$t$','interpreter','latex')
        %title('$\theta$','interpreter','latex')
        title('$local \vec b$','interpreter','latex')
        set(gca,'fontsize',12)
        
        subplot(3,2,4)
        plot(cumsum(theta(1:i)),'linewidth',2)
        xlabel('$t$','interpreter','latex')
        title('cumulative $\theta$','interpreter','latex')
        set(gca,'fontsize',12)
        
        
        subplot(3,2,5)
        scatter(b(1,1:i),b(2,1:i),[],1:i,'linewidth',2)
        xlabel('$b_1(t)$','interpreter','latex')
        ylabel('$b_2(t)$','interpreter','latex')
        title('$\vec b(t)$','interpreter','latex')
        set(gca,'fontsize',12)
        
        
        
        subplot(3,2,6)
        bt(:,i) = R*b(:,i);
        for tt = i-1:-1:1
            R = [cos(theta(tt)), sin(theta(tt)); -sin(theta(tt)), cos(theta(tt))];
            bt(:,i) = R* (bt(:,i) + b(:,tt));
        end
        
       
        
        scatter(bt(1,1:i),bt(2,1:i),[],1:i,'linewidth',2)
        xlabel('$b_1(t)$','interpreter','latex')
        ylabel('$b_2(t)$','interpreter','latex')
        title('cumulative $\vec b(t)$','interpreter','latex')
        set(gca,'fontsize',12)
        set(gcf,'position',[40 154 489 621])
        pause(0.1)
        
        writeVideo(writerObj, getframe(figure(3)));
    end
end



cd ../../../FreezingMethod

if make_movie
    close(writerObj);
end

%%


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



if make_movie
    save([movie_name(1:end-4),'.mat'],'theta','b');
end



