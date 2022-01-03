close all
init = 0.01;
last = 50;
N = 1;
 for i = 1:100
    f = 0.03;
    R = 100;
    %m =250;
    %c =25;
    m = 200 + (300-200).*rand(N,1);
    c = 20 + (30-20).*rand(N,1);
    alpha = 8;
    beta =8;
    beta_d = 8;
    tf =1000;
    eta_IC = [105, -5]';
    eta_dot_IC = [0, 0]';
    theta_hat_IC = [init + (last+init)*rand(N,1), init + (last+init)*rand(N,1)]';
    %theta = [200 + (300-200).*rand(N,1), 20 + (30-20).*rand(N,1) ]';
    theta = [m, c]';
    Gamma = [0.01, 0.01];
    %Gamma = 0.1;
    Gamma = diag(Gamma);
    KCL = 0.01*Gamma;
    dwell_time = 15;
    ts = 0.001;
    t = [0:ts:tf];
    try
        out = sim('finalt3') 
%         list_P=[];
        %error = out.e.Data
    
        f = 0.03;
        R = 100;
        dwell_time = 15;
    %     PX = reshape(out.e_hat.Data(1,:), [], 1);
    %     PPX(i) = PX;
        %PXn = norm(PX);
        %list_P = [list_P,PX];
        e_x(i,:) = out.e.Data(1,:,:);
        e_y(i,:) = out.e.Data(2,:,:);
        r_x(i,:) = out.r.Data(1,:,:);
        r_y(i,:) = out.r.Data(2,:,:);
        es_x(i,:) = out.e_hat.Data(1,:,:);
        es_y(i,:) = out.e_hat.Data(2,:,:);
        es_rx(i,:) = out.r_hat.Data(1,:,:);
        es_ry(i,:) = out.r_hat.Data(2,:,:);
%         er_x(i,:) = [out.e.Data(1,:,:), out.r.Data(1,:,:)];
%         er_y(i,:) = [out.e.Data(2,:,:), out.r.Data(2,:,:)];
%         ehrh_x(i,:) = [out.e_hat.Data(1,:,:), out.r_hat.Data(1,:,:)];
%         ehrh_y(i,:) = [out.e_hat.Data(2,:,:), out.r_hat.Data(2,:,:)];
        Tau_x(i,:) = out.tau.Data(1,:,:);
        Tau_y(i,:) = out.tau.Data(2,:,:);
        ff_x(i,:) = out.ff.Data(1,:,:);
        ff_y(i,:) = out.ff.Data(2,:,:);
        fb_x(i,:) = out.fb.Data(1,:,:);
        fb_y(i,:) = out.fb.Data(2,:,:);
        fbes_x(i,:) = out.fbes.Data(1,:,:);
        fbes_y(i,:) = out.fbes.Data(2,:,:);
        theta_x(i,:) = out.theta_tilde.Data(1,:,:);
        theta_y(i,:) = out.theta_tilde.Data(2,:,:);
%         Px(i,:) = out.e_hat.Data(:,1)';
%         Py(i,:) = out.e_hat.Data(:,2)';
%      
%         Rx(i,:) = out.r.Data(:,1)';
%         Ry(i,:) = out.r.Data(:,2)';
%     
%         EEx(i,:) = out.e.Data(:,1)';
%         EEy(i,:) = out.e.Data(:,2)';
%     
%         REx(i,:) = out.r_hat.Data(:,1)';
%         REy(i,:) = out.r_hat.Data(:,2)';
%     
    end
 end  

    % for jj = 1:10
    %      %PXn(jj) = sum(Px(:,jj).^2)^0.5;
    %      %PYn(jj) = sum(Py(:,jj).^2)^0.5;
    % % %     RXn(jj) = sum(Rx(:,jj).^2)^0.5;
    % % %     EPn(jj) = sum(Ep(:,jj).^2)^0.5;
    % end
     Nex = sqrt(sum(e_x.^2)).^2;
     Ney = sqrt(sum(e_y.^2)).^2;
     Nrx = sqrt(sum(r_x.^2)).^2;
     Nry = sqrt(sum(r_y.^2)).^2;
     Nes_x = sqrt(sum(es_x.^2)).^2;
     Nes_y = sqrt(sum(es_y.^2)).^2;
     Nes_rx = sqrt(sum(es_rx.^2)).^2;
     Nes_ry = sqrt(sum(es_ry.^2)).^2;
     Tau_xx = sqrt(sum(Tau_x.^2)).^2;
     Tau_yy = sqrt(sum(Tau_y.^2)).^2;
     ff_xx = sqrt(sum(ff_x.^2)).^2;
     ff_yy = sqrt(sum(ff_y.^2)).^2;
     fb_xx = sqrt(sum(fb_x.^2)).^2;
     fb_yy = sqrt(sum(fb_y.^2)).^2;
     fbes_xx = sqrt(sum(fbes_x.^2)).^2;
     fbes_yy = sqrt(sum(fbes_y.^2)).^2;
     theta_xx = sqrt(sum(theta_x.^2)).^2;
     theta_yy = sqrt(sum(theta_y.^2)).^2;
%      er_xx = sqrt(sum(er_x.^2)).^2;
%      er_yy = sqrt(sum(er_y.^2)).^2;
%      ehrh_xx = sqrt(sum(ehrh_x.^2)).^2;
%      ehrh_yy = sqrt(sum(ehrh_y.^2)).^2;
     time = linspace(0,10,size(Nex,2));
     figure(1)
     hold on
     plot(time,Nex)
     xlabel('Time', fontsize=12)
     ylabel('Norm of error (x)', fontsize=12)
     plot(time,Ney)
     xlabel('Time', fontsize=12)
     ylabel('Norm of error (Y)', fontsize=12)
     plot(time,Nrx)
     xlabel('Time', fontsize=12)
     ylabel('Norm of filtered Tracking error (x)', fontsize=12)
     plot(time,Nry)
     xlabel('Time', fontsize=12)
     ylabel('Norm of filtered tracking error (y)', fontsize=12)
     xlim([0,5]);
     legend('Error Norm - X', 'Error Norm - Y', 'filtered tracking error - X', 'filtered tracking error - Y')
     grid on
     sgtitle('Norm of error and filtered tracking error for x and y ')

     figure(2)
     hold on
     plot(time,Nes_x)
     xlabel('Time', fontsize=12)
     ylabel('Norm of estimation error (x)', fontsize=12)
     plot(time,Nes_y)
     xlabel('Time', fontsize=12)
     ylabel('Norm of estimation error (y)', fontsize=12)
     plot(time,Nes_rx)
     xlabel('Time', fontsize=12)
     ylabel('Norm of estimated filtered tracking error (x)', fontsize=12)
     plot(time,Nes_ry)
     xlabel('Time', fontsize=12)
     ylabel('Norm of estimated filtered tracking error (y)', fontsize=12)
     xlim([0,5]);
     legend('Estimated Error Norm - X', 'Estimated Error Norm - Y', 'Estimated filtered tracking error - X', 'Estimated filtered tracking error - Y')
     grid on
     sgtitle('Norm of estimated error and filtered tracking error for x and y ')


     figure(3)
     hold on
     plot(time,Tau_xx)
     plot(time,Tau_yy)
     xlabel('Time', fontsize=12)
     ylabel('Tau', fontsize=12)


     %xlim([0,3]);
     %ylim([0,1]);
     grid on
     sgtitle('Norm of Total input ')

     figure(4)
     hold on
     plot(time,ff_xx)
     plot(time,ff_yy)
     xlabel('Time', fontsize=12)
     ylabel('feedforward', fontsize=12)
     legend('feedforward x', 'feedforward y')
     %xlim([0,3]);
     %ylim([0,7e+6]);
     grid on
     sgtitle('Norm of feedforward portion for x and y ')

     figure(5)
     hold on
     plot(time,fb_xx)
     plot(time,fb_yy)
     xlabel('Time', fontsize=12)
     ylabel('feedback', fontsize=12)
     %legend('feedback_x', 'feedback_y')
     legend('feedback for x', 'feedback for y');
     %xlim([0,5]);
     %ylim([0,7e+6]);
     sgtitle('Norm of feedback portion for x and y ')
     grid on

     figure(6)
     hold on
     plot(time,fbes_xx)
     plot(time,fbes_yy)
     xlabel('Time', fontsize=12)
     ylabel('feedback estimation', fontsize=12)
     %legend('feedback_x', 'feedback_y')
     legend('estimated feedback for x','estimated feedback for y');
     %xlim([0,5]);
     %ylim([0,7e+6]);
     sgtitle('Norm of estimated feedback portion for x and y ')
     grid on

     figure(7)
     hold on
     plot(time,theta_xx)
     plot(time,theta_yy)
     xlabel('Time', fontsize=12)
     ylabel('Parametric Errors', fontsize=12)
     legend('Parametric Error - x', 'Parametric Error - y')
    %xlim([0,3]);
    %ylim([0,7e+6]);
     grid on
     sgtitle('Norm of Parametric Error for x and y ')

%     out = sim('finalt3')
%     
%     
%     list_P=[];
% 
%     f = 0.03;
%     R = 100;
%     dwell_time = 15;
%     Px = reshape(out.e_hat.Data(1,:), [], 1);
%     PXn = norm(Px);
%     list_P = [list_P,PXn];
% %    Px(i,:,:) = out.e_hat.Data(1,:,:);
% %    Py(i,:) = out.e_hat.Data(:,2)';
% 
% %    Ex(i,:) = e(:,1)';
% %    Ey(i,:) = e(:,2)';
% % 
% %    Rx(i,:) = r(:,1)';
% %    Ry(i,:) = r(:,2)';
% % 
% %    EEx(i,:) = etahat(:,1)';
% %    EEy(i,:) = etahat(:,2)';
% % 
% %    REx(i,:) = r_hat(:,1)';
% %    REy(i,:) = r_hat(:,2)';
% 
% % end
% % for jj = 1:size(Px,2)
% %     PXn(jj) = sum(Px(:,jj).^2)^0.5;
% %     PYn(jj) = sum(Py(:,jj).^2)^0.5;
% % %     RXn(jj) = sum(Rx(:,jj).^2)^0.5;
% % %     EPn(jj) = sum(Ep(:,jj).^2)^0.5;
% % NPx = sqrt(sum(Px.^2)).^2;
% % NPy = sqrt(sum(Py.^2)).^2;
% % time = out.e_hat.Time(:,1)';
% % Npoint = 1000;
figure(8)
subplot(2,1,1)
    hold on
    plot(out.e.Data(1,:))
    plot(out.e.Data(2,:))
    xlabel('Time',fontsize =12);
    ylabel('e',fontsize =12);
    grid on
    subplot(2,1,2)
    hold on
    plot(out.r.Data(1,:))
    plot(out.r.Data(2,:))
    %plot(out.e.Time(2,:), out.e.Data(2,:))
    xlabel('Time',fontsize =12);
    ylabel('r',fontsize =12);
    grid on
    sgtitle('Error and Filtered Tracking Error')
figure(9)
    subplot(1,3,1)
    hold on
    plot(out.eta.Data(1,:), out.eta.Data(2,:),'b')
    plot(out.eta_d.Data(1,:), out.eta_d.Data(2,:),'--r')
    xlabel('x position',fontsize =12);
    ylabel('y position',fontsize =12);
    legend('state','Desired state');
    % xlim([-50,50]);
    % ylim([-50,50]);
    grid on
    subplot(1,3,2)
    hold on
    plot(out.eta_dot.Data(1,:), out.eta_dot.Data(2,:),'b')
    plot(out.eta_d_dot.Data(1,:), out.eta_d_dot.Data(2,:),'--r')
    xlabel('x velocity',fontsize =12);
    ylabel('y velocity',fontsize =12);
    % xlim([-50,50]);
    % ylim([-50,50]);
    legend('state','Desired state');
    grid on
    subplot(1,3,3)
    hold on
    plot(out.eta_dotdot.Data(1,:), out.eta_dotdot.Data(2,:),'b')
    plot(out.eta_d_dotdot.Data(1,:), out.eta_d_dotdot.Data(2,:),'--r')
    xlabel('x acceleration',fontsize =12);
    ylabel('y acceleration',fontsize =12);
    legend('state','Desired state');
    % xlim([-50,50]);
    % ylim([-50,50]);
    grid on
%     sgtitle('Tracjectory when feedback does not exist all the time and \n parameters are unknown')
    sgtitle({['Tracjectory when feedback does not exist all the time and'] ['parameters are known']},fontsize =13)
%     figure(2)
%     subplot(2,1,1)
%     hold on
%     plot(out.e_hat.Data(1,:))
%     plot(out.e_hat.Data(2,:))
%     xlabel('Time',fontsize =12);
%     ylabel('e_hat',fontsize =12);
%     grid on
%     subplot(2,1,2)
%     hold on
%     plot(out.r_hat.Data(1,:))
%     plot(out.r_hat.Data(2,:))
%     %plot(out.e.Time(2,:), out.e.Data(2,:))
%     xlabel('Time',fontsize =12);
%     ylabel('r_hat',fontsize =12);
%     grid on
%     sgtitle('Estimation Error and Filtered Tracking Error')
% 
% % figure(3)
% % subplot(211)
% % hold on
% % plot(list_P)
% % % plot(PXn,'--')
% % % plot(PYn,'o')
% % %plot(NPx(1:Npoint:end),'or','Linewidth',3);
% % xlabel('Time');
% % ylabel('Position Error');
% % %ylim([0,0.1])
% % grid on
% 


% figure(4)
% hold on
% plot(out.theta_tilde.Data(1,:))
% plot(out.theta_tilde.Data(2,:))
% xlabel('Time',fontsize =12);
% ylabel('theta tilde',fontsize =12);
% grid on
% sgtitle('Parametric Error')
% 
% figure(5)
% subplot(2,2,1)
% hold on
% plot(out.tau.Data(1,:))
% plot(out.tau.Data(2,:))
% xlabel('Time',fontsize =12);
% ylabel('Total Input',fontsize =12);
% grid on
% subplot(2,2,2)
% hold on
% plot(out.ff.Data(1,:))
% plot(out.ff.Data(2,:))
% xlabel('Time',fontsize =12);
% ylabel('Feedforward',fontsize =12);
% grid on
% subplot(2,2,3)
% hold on
% plot(out.fb.Data(1,:))
% plot(out.fb.Data(2,:))
% xlabel('Time',fontsize =12);
% ylabel('feedback',fontsize =12);
% grid on
% subplot(2,2,4)
% hold on
% plot(out.fbes.Data(1,:))
% plot(out.fbes.Data(2,:))
% xlabel('Time',fontsize =12);
% ylabel('Feedback Estimation',fontsize =12);
% grid on

