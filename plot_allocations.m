%%
clc;
clear all;


Wa =        7; % width for allocation subplot
Wr =        1; % width for realibility subplot
Ws =        2; % width for space

W_v =       [Wa,Ws,repmat([Wa Wr Ws],1,4)];
W_cumsum =  cumsum(W_v);
Wt =        sum(W_v); % total width
W_start =   [1,1+W_cumsum(1:end-1)];
W_end =     W_cumsum;

mis0adapt0 =  load(['online_cp_over_iters_mis0_adapt0.mat']);
mis0adapt1 =  load(['online_cp_over_iters_mis0_adapt1.mat']);
mis1adapt0 =  load(['online_cp_over_iters_mis1_adapt0.mat']);
mis1adapt1 =  load(['online_cp_over_iters_mis1_adapt1.mat']);
    
assert(all(mis0adapt0.G_alc_f(:)==mis0adapt1.G_alc_f(:))); % making sure all have the same traffic
assert(all(mis0adapt0.G_alc_f(:)==mis1adapt0.G_alc_f(:)));
assert(all(mis0adapt0.G_alc_f(:)==mis1adapt1.G_alc_f(:)));

N = mis0adapt0.num_frames;
YLim = [N-200,N];

figure;
sp1 = subplot(1,Wt, W_start( 1) : W_end( 1) ); imagesc(mis0adapt0.G_alc_f); set(gca,'YLim',YLim);                                 title('URLLC generation','interpreter','latex','FontSize',18); xlabel('Slots','interpreter','latex','FontSize',18); ylabel('Frames','interpreter','latex','FontSize',18);
sp3 = subplot(1,Wt, W_start( 3) : W_end( 3) ); imagesc(mis0adapt0.U_alc_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); title('Conventional scheduler','interpreter','latex','FontSize',18);
sp4 = subplot(1,Wt, W_start( 4) : W_end( 4) ); imagesc(mis0adapt0.r_mux_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); 
sp6 = subplot(1,Wt, W_start( 6) : W_end( 6) ); imagesc(mis0adapt1.U_alc_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); title('CP-based scheduler'    ,'interpreter','latex','FontSize',18);
sp7 = subplot(1,Wt, W_start( 7) : W_end( 7) ); imagesc(mis0adapt1.r_mux_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); 
sp9 = subplot(1,Wt, W_start( 9) : W_end( 9) ); imagesc(mis1adapt0.U_alc_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); title('Conventional scheduler','interpreter','latex','FontSize',18);
sp10= subplot(1,Wt, W_start(10) : W_end(10) ); imagesc(mis1adapt0.r_mux_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]);
sp12= subplot(1,Wt, W_start(12) : W_end(12) ); imagesc(mis1adapt1.U_alc_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); title('CP-based scheduler'    ,'interpreter','latex','FontSize',18);
sp13= subplot(1,Wt, W_start(13) : W_end(13) ); imagesc(mis1adapt1.r_mux_f); set(gca,'YLim',YLim); set(gca,'xtick',[],'ytick',[]); title('reliability indicator' ,'interpreter','latex','FontSize',18);

colormap(sp1, cool);
colormap(sp3, cool);
colormap(sp6, cool);
colormap(sp9, cool);
colormap(sp12,cool);

% some printings
covrg_mis0 = [mis0adapt0.covrg_time_avg(end)      ; mis0adapt1.covrg_time_avg(end)]
effic_mis0 = [mis0adapt0.effic_eMBB_time_avg(end) ; mis0adapt1.effic_eMBB_time_avg(end)]
covrg_mis1 = [mis1adapt0.covrg_time_avg(end)      ; mis1adapt1.covrg_time_avg(end)]
effic_mis1 = [mis1adapt0.effic_eMBB_time_avg(end) ; mis1adapt1.effic_eMBB_time_avg(end)]

set(gcf,'Position',[0 0 1240 680]);