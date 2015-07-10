%% this code is a wrapper to perform simulations with a BinaryTree model
% the signal strength is normalized so that the oracle achieves 1% error
cd('C:\Git/high-dim-risk-experiments/Experiments/Ridge/BinaryTree');
addpath '../../../Code/'

%%  normalize the signal strength to have the same oracle mis-classification error

gamma = [0.5; 1; 2];
alpha = [0.5; 1; 2];
depth = [2; 4; 6]; %p = 2^depth - 1

%% loop over parameters
a = {'-','--',':','-.'};
n_lambda = 50;
for d = 1:length(depth)
    for j=1:length(gamma)
        for k=1:length(alpha)
            [lambda,estim_risk,pred_risk,lambda_th,estim_risk_th,pred_risk_th] = ridge_hier_sim(gamma(j),alpha(k),depth(d),n_lambda);
    
            %save results
            figure, hold on
            normalize =  @(lambda)  sqrt(lambda);
            h = plot(normalize(lambda),estim_risk);
            set(h,'Linewidth',3,'LineStyle',a{1})
            h = plot(normalize(lambda_th),estim_risk_th);
            set(h,'Linewidth',3,'LineStyle',a{2})
            xlabel('$\sqrt{\lambda}$','Interpreter','LaTex')
            xlim([min(normalize(lambda)) max(normalize(lambda))]);
            ylim([0,1.1*max(max(estim_risk),max(estim_risk_th))])
            %legend('Empirical Est','Theoretical Est','Location','Best');
            set(gca,'fontsize',20)
            %titlestr = sprintf( 'depth=%.1f gamma=%.1f alpha=%.1f',depth(d),gamma(j),alpha(k));
            %title(titlestr);
            filename = sprintf( './hier_depth_%d_gamma_%.1f_alpha_%.2f_estim.pdf',depth(d),gamma(j),alpha(k));
            save2pdf(filename)
            %saveas(gcf, filename,'pdf');
            fprintf(['Saved Results to ' filename '\n']);
            
            figure, hold on
            h = plot(normalize(lambda),pred_risk);
            set(h,'Linewidth',3,'LineStyle',a{1})
            h = plot(normalize(lambda_th),pred_risk_th);
            set(h,'Linewidth',3,'LineStyle',a{2})
            xlabel('$\sqrt{\lambda}$','Interpreter','LaTex')
            xlim([min(normalize(lambda)) max(normalize(lambda))]);
            ylim([0,1.1*max(max(pred_risk),max(pred_risk_th))])
            %legend('Empirical Pred','Theoretical Pred','Location','Best');
            set(gca,'fontsize',20)
            %titlestr = sprintf( 'depth=%.1f gamma=%.1f alpha=%.1f',depth(d),gamma(j),alpha(k));
            %title(titlestr);
            filename = sprintf( './hier_depth_%d_gamma_%.1f_alpha_%.2f_pred.pdf',depth(d),gamma(j),alpha(k));
            save2pdf(filename)
            %saveas(gcf, filename,'pdf');
            fprintf(['Saved Results to ' filename '\n']);
      
        end
    end
end