%% this code is a wrapper to perform simulations with an exponential model
% the signal strength is normalized so that the oracle achieves 1% error
cd('C:\Git\high-dim-risk-experiments\Experiments\RDA\exponential');
addpath '..\..\..\Code\'

%%  normalize the signal strength to have the same oracle mis-classification error
%run AR(1) simulations for a grid of parameters

%define grids
%exp_lambda = [1/10; 1/4; 1; 4];
%gamma = [0.5; 1; 2];
%oracle_error_rate =0.01; %[0.01 0.02 0.05];
%% test case
gamma = [0.5; 1; 2];
exp_lambda =1;
oracle_error_rate = 0.01;

%% loop over parameters
a = {'-','--',':','-.'};
n_lambda = 100;
n = 5*1e2;
for i=1:length(exp_lambda)
    for j=1:length(gamma)
        for k=1:length(oracle_error_rate)
            %perform experiment
            
            [lambda,risk,lambda_th,risk_th] = run_exp_sim_norm_sig(exp_lambda(i),gamma(j),oracle_error_rate(k),n_lambda,n);
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
            xlabel('$\sqrt{\lambda}$','Interpreter','LaTex')
            xlim([min(normalize(lambda)) max(normalize(lambda))]);
            %set(gca,'XTick',[0:5])
            %set(gca,'XTickLabel',{'0','1','2','3','4','5'})
            ylim([0, 0.5]);
            ylabel('test error');
            %set(gca,'YTick',0.1*[0:5])
            %set(gca,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'})
            
            %legend('Empirical P(err)','Theoretical P(err)','Oracle','Location','Best');
            set(gca,'fontsize',20)
            %titlestr = sprintf( 'l=%.2f gamma=%.1f err=%.2f',exp_lambda(i),gamma(j),oracle_error_rate(k));
            %title(titlestr);
            filename = sprintf( './exp_l_%.2f_gamma_%.1f_err_%.2f.pdf',exp_lambda(i),gamma(j),oracle_error_rate(k));
            save2pdf(filename)
            %saveas(gcf, filename,'pdf');
            fprintf(['Saved Results to ' filename '\n']);
            
        end
    end
end
