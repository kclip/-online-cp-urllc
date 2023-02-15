clear all;
clc;
close all


for b_mismatch_case = 0:1
    for b_is_adaptive = 0:1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% external loops
        if b_mismatch_case==0
            p_plus_pop =        0.16;
            p_minus_pop =       0.16;
            p_plus_hat =        0.02;
            p_minus_hat =       0.02;
            G_num_min_pop =     0;
            G_num_max_pop =     6;
            G_num_min_hat =     0;
            G_num_max_hat =     6;
            lr_theta =          0.1 * b_is_adaptive;
        else
            p_plus_pop =        0.16;
            p_minus_pop =       0.16;
            p_plus_hat =        0.40;
            p_minus_hat =       0.40;
            G_num_min_pop =     0;
            G_num_max_pop =     6;
            G_num_min_hat =     0;
            G_num_max_hat =     6;
            lr_theta =          0.1 * b_is_adaptive;
        end
        
        theta_M =           1;
        alpha =             0.1;
        num_frames =        2000;
        G_num_0 =           2;
        S =                 12;
        L =                 1;
        b_0ind_1rate =      0;
        b_ShowInnerPlots =  0;

        varphi =           @(theta) 0.5 * ( 1 + sin (pi*(max(0,min(1,theta/theta_M))-0.5)) ); % stretching func
        varphi_inv =       @(phi)   theta_M*(0.5+(asin(2*phi-1)/pi)); % inverse to the stretching function within transition region


        rng_seed =          654;
        rng(rng_seed);
        rand_0_to_1_f =     rand(   num_frames,1);
        W_pop_f =           nan(    num_frames,1);
        G_num_f =           nan(    num_frames,1);
        U_alc_f =           nan(    num_frames,S);
        Gamma_int_f =       cell(   num_frames,1);
        G_alc_f =           nan(    num_frames,S);
        miscovered_alc_f =  nan(    num_frames,S);
        r_ind_f =           nan(    num_frames,1);
        r_rate_f =          nan(    num_frames,1);
        r_exact_f =         nan(    num_frames,1);
        r_mux_f =           nan(    num_frames,1);
        theta_f =           nan(    num_frames,1);
        q_f =               zeros(  num_frames,2^S);
        G_alc_f(1,:) =                              0; % G_num_0 slots will be overiden with 1
        G_alc_f(1,find(randperm(S)<=G_num_0)) =     1; % allocating in random G_num_0 slots within the first frame
        U_alc_f(1,:) =                              0; % G_num_max slots will be overiden with 1
        miscovered_alc_f(1,:) =                     0;
        G_num_f(1) =                                G_num_0;
        theta_f(1:2) =                              varphi_inv(alpha); % inverse to the stretching function within transition region
        r_ind_f(1) =                                0;
        r_rate_f(1) =                               0;
        r_exact_f(1) =                              0;
        r_mux_f(1) =                                0;

        
        for ff=2:num_frames      % building the traffic upfront, so we can reproduce using the same seed without the CP randomness of the tie breakers changing the traffic
            % build G_alc_ff
            W_pop_f(ff) =                        (-1)*(rand_0_to_1_f(ff) <     p_minus_pop*(G_num_f(ff-1)>G_num_min_pop))...    % -1 w.p. p^-(G_fm1)
                                                +(+1)*(rand_0_to_1_f(ff) > 1 - p_plus_pop *(G_num_f(ff-1)<G_num_max_pop) )...   % +1 w.p. p^+(G_fm1)
                                                + 0;                                                                            %  0 w.p. 1 - p^+ - p^-
            G_num_f(ff) =                       G_num_f(ff-1) + W_pop_f(ff);
            G_alc_f(ff,:) =                     G_alc_f(ff-1,:); % to be overriden on birth or death events
            switch W_pop_f(ff)
                case 1  % birth
                    unoccupied_vec =            set_alc2vec(1-G_alc_f(ff-1,:));
                    unoccupied_vec =            unoccupied_vec{1};
                    slot_birth =                unoccupied_vec(randi(length(unoccupied_vec)));
                    G_alc_f(ff,slot_birth) =    1;
                case -1 % death
                    occupied_vec =              set_alc2vec(  G_alc_f(ff-1,:));
                    occupied_vec =              occupied_vec{1};
                    slot_death =                occupied_vec(  randi(length(occupied_vec  )));
                    G_alc_f(ff,slot_death) =    0;
            end        
        end

        for ff=2:num_frames
            % build q_f:
            G_num_fminus1 =                     sum(G_alc_f(ff-1,:),2);
            [neighbor_p1_int,...
             neighbor_m1_int,...
             neighbor_00_int] =                 neighbors_vec(      G_alc_f(ff-1,:), G_num_min_hat, G_num_max_hat );
            q_f(ff,1+neighbor_p1_int) =         p_plus_hat *(G_num_fminus1<G_num_max_hat) / length(neighbor_p1_int) ;      % birth      w.p. p_plus_hat
            q_f(ff,1+neighbor_m1_int) =         p_minus_hat*(G_num_fminus1>G_num_min_hat) / length(neighbor_m1_int) ;      % death      w.p. p_minus_hat
            q_f(ff,1+neighbor_00_int) =         (1 - p_plus_hat  * (G_num_fminus1<G_num_max_hat) ...
                                                   - p_minus_hat * (G_num_fminus1>G_num_min_hat)      ) / length(neighbor_00_int)   ;  % stationary; only one set, w.p. 1-p_plus_hat-p_minus_hat
            % naive set predictor Gamma_f:
            %theta_th =                          max(0.0,min(1.0,theta_f(ff)));
            var_phi =                           varphi(theta_f(ff)); % stretching func
            sorted_idx =                        sort_random_tiebreakers(q_f(ff,:),'descend');
            cdf_vec =                           cumsum(q_f(ff,sorted_idx));
            if var_phi==0.0
                U_alc_f(ff,:) =                 ones(1,S); % technical condition. if theta hits 1, we must have the entire set within
            else
                pred_set_len =                  min(find(cdf_vec>=1-var_phi));
                if pred_set_len>1 % randomness can take this by 1, so making sure only if it has 2
                    p_above =                   cdf_vec(pred_set_len);
                    p_below =                   cdf_vec(pred_set_len-1);
                else
                    p_above =                   cdf_vec(pred_set_len);
                    p_below =                   0;
                end
                rand_cutoff =                   (p_below + (p_above-p_below) * rand(1)) < (1-var_phi); % 0 or 1
                pred_set_len =                  pred_set_len-1  + rand_cutoff;
                if pred_set_len>=1
                    Gamma_int_f{ff} =           sorted_idx(1:pred_set_len)-1;  % index back to being a allocation in integer representation
                else
                    Gamma_int_f{ff} =           []; % predict the empty set
                end
                % greedy allcate U_f to serve all G_f in Gamma_f:
                U_alc_f(ff,:) =                 greedy_allocation(Gamma_int_f{ff},L,S);
            end

            % retrieve reliability of allcoation
            [miscovered_alc_f(ff,:), ...
                r_ind_f(ff), ...
                r_rate_f(ff),...
                r_exact_f(ff)] =                coverage_reliability (G_alc_f(ff,:) , U_alc_f(ff,:) , L);
            % update calibration parameter
            r_mux_f(ff) =                       (  b_0ind_1rate) * r_rate_f(ff) +...
                                                (1-b_0ind_1rate) * r_ind_f( ff) ;
            theta_f(ff+1) =                     theta_f(ff) + lr_theta*(   r_mux_f(ff) -(1-alpha)) ;
        end

        effic_eMBB_f =                          (S -sum(U_alc_f,2))/S;
        effic_eMBB_time_avg =                   cumsum(effic_eMBB_f) ./ (1:num_frames).';
        covrg_time_avg =                        cumsum(r_mux_f     ) ./ (1:num_frames).';

        %%

        if b_ShowInnerPlots
            figure; plot(effic_eMBB_f,'-r'); hold on;
            plot(effic_eMBB_time_avg,'b','LineWidth',2); title('eMBB efficiency');

            figure; plot(G_num_f); title('G num'); axis([1,num_frames,G_num_min_pop-1,G_num_max_pop+1]);
            figure; plot(covrg_time_avg); title('coverage'); hold on;
            plot(r_mux_f,'rx');
            plot([1,num_frames],(1-alpha)*[1,1],'--g'); set(gca,'YLim',[0,1]);

            figure; 
            subplot(1,4,1); imagesc(G_alc_f);           axis([1,S,num_frames-100,num_frames]); title('Generated packets'); xlabel('Slots within frame'); ylabel('frames');
            subplot(1,4,2); imagesc(U_alc_f);           axis([1,S,num_frames-100,num_frames]); title('Allocation');
            subplot(1,4,3); imagesc(U_alc_f+3*G_alc_f); axis([1,S,num_frames-100,num_frames]); title('U on top of G');  colormap cool;
            subplot(1,4,4); imagesc(miscovered_alc_f);  axis([1,S,num_frames-100,num_frames]); title('L-miscovered'); %colorbar;

            figure; plot(theta_f); title('theta');
        end

        save(['online_cp_over_iters_mis',num2str(b_mismatch_case),'_adapt',num2str(b_is_adaptive)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% external loops
    end
end

plot_allocations

%%%%%%%%%%%%%% AUX functions %%%%%%%%%%%%%% 

function [neighbor_p1_int,neighbor_m1_int,neighbor_00_int] = neighbors_vec(G_alc_fminus1, G_num_min, G_num_max)
    G_num_fminus1 =                     sum(G_alc_fminus1);
    S =                                 size(G_alc_fminus1,2);
    eye_S =                             eye(S);
    % birth: total of S-G_num_fminus1 subsets
    if (G_num_fminus1<G_num_max)
        neighbor_p1_alc =               G_alc_fminus1 + (eye_S(:,find(G_alc_fminus1==0))).'; % bsx row+mat
        neighbor_p1_int =               set_vec2int(set_alc2vec(neighbor_p1_alc));
    else
        neighbor_p1_int =               [];
    end
    % death: total of   G_num_fminus1 subsets. 
    if (G_num_fminus1>G_num_min)
        neighbor_m1_alc =               G_alc_fminus1 - (eye_S(:,find(G_alc_fminus1==1))).'; % bsx row-mat
        neighbor_m1_int =               set_vec2int(set_alc2vec(neighbor_m1_alc));
    else
        neighbor_m1_int =               [];
    end
    % stationary
    neighbor_00_int =                   set_vec2int(set_alc2vec(G_alc_fminus1));        % return only 1 subset, the same as the input 
end

function set_int = set_vec2int(set_vec)% works over cell array of vectors
    set_int =                   zeros(length(set_vec),1);    
    for ii=1:length(set_vec)
        if isempty(set_vec{ii})
            set_int(ii) =       0;
        else
            set_int(ii) =       sum(2.^(set_vec{ii}-1));
        end
    end
end

function set_alc = set_int2alc(set_int,S) % works over a column vector (or a scalar)
    set_alc =           fliplr(dec2bin(set_int,S))=='1';
end

function set_vec = set_alc2vec(set_alc) % works over matrix (or a row vector)
    set_vec =           cell(size(set_alc,1),1);
    for ii=1:size(set_alc,1)
        set_vec{ii} =   find(set_alc(ii,:));  % must act on rows and not on matrix...
    end
end

function U_alc = greedy_allocation(Gamma_int,L,S)
    U_alc =                                         zeros(1,S);    
    if ~isempty(Gamma_int)
        Gamma_alc =                                 set_int2alc(Gamma_int,S); % copy that is changes later
        for s=S:-1:1 % backwards
            if any(Gamma_alc(:,s))
                U_alc(s) =                          1;  % add s
                set_vec =                           set_alc2vec(Gamma_alc); % reevaluate since it may change in the iteration
                for gg=1:size(Gamma_alc,1) % go over all sets within Gamma
                    covered_slot =                  max(intersect(s-(L:-1:0),set_vec{gg}));
                    Gamma_alc(gg,covered_slot) =    0; % remove the latest generated slot that is L-covered
                end
            end
        end
    end
end

function [miscovered_alc, r_ind , r_rate, r_exact] = coverage_reliability (G_alc , U_alc , L)
    G_vec_orig =                        set_alc2vec(G_alc);
    G_vec_orig =                        G_vec_orig{1};
    G_vec =                             G_vec_orig;
    G_num =                             length(G_vec);
    U_vec =                             set_alc2vec(U_alc);
    U_vec =                             U_vec{1};
    for s=fliplr(U_vec)
        covered_slot =                  max(intersect(G_vec,s-(L:-1:0)));
        G_vec =                         setdiff(G_vec,covered_slot); % "serve", and remove it
    end
    if G_num>0
        r_rate =                        1 - length(G_vec)/G_num; % ratio of L-covered
        r_ind =                         (r_rate==1); % boolean = if all generated traffic was served
    else
        r_rate =                        1; % if no packets were generated, we treat as if full L-coverage
        r_ind =                         1; % boolean = if all generated traffic was served
    end
    miscovered_alc =                    set_int2alc(set_vec2int({G_vec}),size(G_alc,2));
    r_exact =                           any(all(U_alc==G_alc)); % is U_alc a row in G_alc
end

function sorted_idx = sort_random_tiebreakers(x,direction) % x must be a row vector
    assert(             isrow(x));
    [v,idx] =           sort(x,direction);
    idxc =              cumsum([1 logical(diff(v))]);
    sorted_idx =        cell2mat(accumarray(idxc',idx',[],@(x){x(randperm(numel(x)))}));
end

