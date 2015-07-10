%% this code is a wrapper to perform simulations with an exponential model
% the signal strength is normalized so that the oracle achieves 1% error
cd('C:\Git/high-dim-risk-experiments/Experiments/RDA/hierarch');
addpath '../../../Code/'

%%  normalize the signal strength to have the same oracle mis-classification error
%run AR(1) simulations for a grid of parameters


gamma = [0.5; 1; 2];
oracle_error_rate = [0.01; 0.05];
depth = [9; 10; 11]; %p = 2^depth - 1

%% loop over parameters
a = {'-','--',':','-.'};
n_lambda = 100;
for d = 1:length(depth)
    for j=1:length(gamma)
        for k=1:length(oracle_error_rate)
            %perform experiment
            
            [lambda,risk,lambda_th,risk_th] = run_hier_sim_norm_sig(gamma(j),oracle_error_rate(k),n_lambda,depth(d));
            %%
            
            %save results
            figure, hold on
            %normalize =  @(lambda)  (2*((lambda.*(lambda+4)).^(1/4))./(lambda.^(1/2)+(lambda+4).^(1/2))).^2;
            normalize =  @(lambda)  sqrt(lambda);
            h = plot(normalize(lambda),risk);
            set(h,'Linewidth',3,'LineStyle',a{1})
            h = plot(normalize(lambda_th),risk_th);
            set(h,'Linewidth',3,'LineStyle',a{2})
            h = plot(normalize(lambda_th),oracle_error_rate(k)*ones(length(lambda_th),1));
            set(h,'Linewidth',3,'LineStyle',a{3})
            xlabel('$$\sqrt{\lambda}$$','Interpreter','latex','FontSize',20);
            xlim([min(normalize(lambda)) max(normalize(lambda))]);
            ylim([0, 0.5]);
            ylabel('test error');
            %legend('Empirical P(err)','Theoretical P(err)','Oracle','Location','Best');
            set(gca,'fontsize',20)
            %titlestr = sprintf( 'l=%.2f gamma=%.1f err=%.2f',exp_lambda(i),gamma(j),oracle_error_rate(k));
            %title(titlestr);
            filename = sprintf( './hier_depth_%d_gamma_%.1f_err_%.2f.pdf',depth(d),gamma(j),oracle_error_rate(k));
            save2pdf(filename)
            %saveas(gcf, filename,'pdf');
            fprintf(['Saved Results to ' filename '\n']);
            
        end
    end
end