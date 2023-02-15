clear all;
clc;
close all;

tic

alpha =                                 0.01;
lr_theta =                              0.05; % learning rate for updating theta
v_p_pop =                               0.08*2.^(-5:0.5:0);
v_p_hat =                               0.32*2.^(-7:0);
v_b_0ind_1rate =                        [0,1];
[m_p_pop,m_p_hat,m_b_0ind_1rate] =      meshgrid(v_p_pop, v_p_hat, v_b_0ind_1rate);
rng_seed =                              123;
A =                                     ones(size(m_p_pop)); % for making scalars being the same matrix dimensions

[m_effic_adapt, m_covrg_adapt] =        arrayfun(...
                                            @cp_urllc_one_trial ,...
                                            lr_theta        * A ,...
                                            m_p_pop             ,...
                                            m_p_pop             ,...
                                            m_p_hat             ,...
                                            m_p_hat             ,...
                                            m_b_0ind_1rate      ,...
                                            alpha           * A ,...
                                            rng_seed        * A);
%%
[m_effic_nonadapt, m_covrg_nonadapt] =  arrayfun(...
                                            @cp_urllc_one_trial ,...
                                            lr_theta    * 0 * A ,...
                                            m_p_pop             ,...
                                            m_p_pop             ,...
                                            m_p_hat             ,...
                                            m_p_hat             ,...
                                            m_b_0ind_1rate      ,...
                                            alpha           * A ,...
                                            rng_seed        * A);
toc
save main_online_sweep_pop_and_hat
%%

if ~exist('m_effic_adapt','var')
    load('main_online_sweep_pop_and_hat.mat');
end

%close all;
effic_levels_0_color =                      0.000:0.001:1.000;
effic_levels_0_text =                       0.1:0.1:1.0;
covrg_levels_0_dashed =                     1-alpha;
covrg_levels_0_color =                      0.90:0.001:1.00;
covrg_levels_0_text =                       setdiff([0.90:0.01:0.98,0.985:0.005:1.00],covrg_levels_0_dashed);
effic_levels_1_color =                      effic_levels_0_color;
effic_levels_1_text =                       effic_levels_0_text;
covrg_levels_1_dashed =                     1-alpha;
covrg_levels_1_color =                      covrg_levels_0_color;
covrg_levels_1_text =                       setdiff(0.97:0.001:1.00,covrg_levels_1_dashed);


res =                                       0.002;

m_effic_adapt =                             res * round(1/res * m_effic_adapt);
m_covrg_adapt =                             res * round(1/res * m_covrg_adapt);
m_effic_nonadapt =                          res * round(1/res * m_effic_nonadapt);
m_covrg_nonadapt =                          res * round(1/res * m_covrg_nonadapt);

for ii=1:length(v_b_0ind_1rate)
    figure;
    dashed = contour(m_p_pop(:,:,ii), m_p_hat(:,:,ii), m_covrg_nonadapt(:,:,ii),covrg_levels_0_dashed*[1,1]);
    close;
    figure; axis tight; axis equal square;
    sp1 = subplot(2,2,1); plot_contour(m_p_pop(:,:,ii), m_p_hat(:,:,ii), m_covrg_nonadapt(:,:,ii),covrg_levels_0_color,covrg_levels_0_text,covrg_levels_0_dashed ); title('URLLC reliability rate' ,  'Interpreter','latex','FontSize',18);  colorbar;
    sp2 = subplot(2,2,2); plot_contour(m_p_pop(:,:,ii), m_p_hat(:,:,ii), m_effic_nonadapt(:,:,ii),effic_levels_0_color,effic_levels_0_text,[]                    ); title('Empirical eMBB efficiency','Interpreter','latex','FontSize',18'); colorbar; plot(dashed(1,2:end)+1j*dashed(2,2:end),'r--','LineWidth',2);
    sp3 = subplot(2,2,3); plot_contour(m_p_pop(:,:,ii), m_p_hat(:,:,ii), m_covrg_adapt(   :,:,ii),covrg_levels_1_color,covrg_levels_1_text,[] );                    title('URLLC reliability rate'  , 'Interpreter','latex','FontSize',18);  colorbar;
    sp4 = subplot(2,2,4); plot_contour(m_p_pop(:,:,ii), m_p_hat(:,:,ii), m_effic_adapt(   :,:,ii),effic_levels_1_color,effic_levels_1_text,[]);                     title('Empirical eMBB efficiency','Interpreter','latex','FontSize',18);  colorbar;
    xlabel('Ground-truth parameter $p$','Interpreter','latex','FontSize',16); ylabel('Predicotr parameter $\hat{p}$','Interpreter','latex','FontSize',16);
    colormap(sp1, spring);
    colormap(sp3, spring);
    set(gcf,'Position',[0 0 960 600]);
end

function plot_contour(X,Y,Z,levels_color,levels_text,levels_dashed)
    contourf(X,Y,Z,levels_color,'LineColor','none');%,'ShowText','off');
    hold on;
    contour(X,Y,Z,levels_text,'LineColor','k','ShowText','on');
    if ~isempty(levels_dashed)
        if isscalar(levels_dashed)
            levels_dashed =     levels_dashed*[1,1]; % for visualization purposes
        end
        contour(X,Y,Z,[levels_dashed],'r--','LineWidth',2,'ShowText','on');
    end
    set(gca,'CLim',[levels_color(1),levels_color(end)])
    xlabel('$p$','Interpreter','latex','FontSize',16); ylabel('$\hat{p}$','Interpreter','latex','FontSize',16);
    axis tight; set(gca,'XScale','log');set(gca,'YScale','log');
    xticks(X(1,1:2:end)); xticklabels(num2str(X(1,1:2:end).'));
    yticks(Y(1:2:end,1)); yticklabels(num2str(Y(1:2:end,1)  ));
end


