function cfg = config()
% config.m  --  Central configuration for the frontoparietal GAMMA connectivity
% pilot (Pilot 2). Dataset: OpenNeuro ds003702 "Social Memory Cuing"
% (Gregory, Wang & Kessler; VR social/non-social cueing, visual WM task,
% 64-ch eego sports system).
%
% VALUES BELOW ARE NOW SET from: (a) s00_inventory.m run on sub-01,
% (b) ReadMe-EventCodes.txt (authors' trigger documentation),
% (c) the authors' trialfun_TakesInputsCueToTendRespTs.m and
%     TrialInformationEEGProcessedData.txt (trialinfo column definitions),
% (d) direct inspection of a real trialinfo matrix (224x17, sub-01).
%
% KEY ARCHITECTURAL DECISION (different from Pilot 1): trial definition
% (lock samples, cue type, validity, location-match, accuracy) is NOT
% re-derived from raw .vmrk marker pairs in our own code. Instead we load
% the authors' own `trialinfo` matrix directly from each subject's
% sub-XX\ProcessedData\data_ica.mat (only trialinfo -- NOT the 0.5-36 Hz
% filtered EEG data itself, which is unusable for gamma). trialinfo
% columns 8-15 are ABSOLUTE sample indices
% into the same raw continuous file we load ourselves via pop_loadbv, so
% they line up directly with OUR independently-preprocessed (wideband)
% continuous data, AS LONG AS we do not resample (see cfg.resample_hz).
%
% Run order:  s00_inventory -> s01_preprocess -> s02_epoch (trialinfo-based)
%             -> s03_connectivity -> s04_stats -> s05_figures

% ---------------- PATHS (EDIT to your machine) ----------------
cfg.bids_dir    = fullfile('..','data','ds003702');                          % EDIT to your local path
cfg.deriv_dir   = fullfile('..','derivatives','fpg');                        % EDIT to your local path (fpg = FrontoParietal Gamma)
% CONFIRMED: unlike the ReadMeDerivatives.txt description (which implied a
% separate derivatives_ProcessedData/EEGDataTableStudy/sub-XX/... tree),
% on this machine data_ica.mat actually lives directly under each subject's
% own folder inside cfg.bids_dir: sub-XX\ProcessedData\data_ica.mat.
% cfg.authors_deriv_file(subj) below builds that path per subject.
cfg.authors_deriv_file = @(subj) fullfile(cfg.bids_dir, subj, 'ProcessedData', 'data_ica.mat');
cfg.eeglab_dir  = fullfile('..','..','eeglab2023.0');   % EDIT to your EEGLAB install

% ---------------- DATASET SUBSET (for a fast first pass) ----------------
cfg.subset_n = [];   % [] = process all subjects for the final run

% ---------------- PREPROCESSING ----------------
cfg.resample_hz = 500;   % CONFIRMED native rate (s00_inventory). Deliberately
                          % NOT resampling, to keep trialinfo's absolute sample
                          % indices valid without rescaling. If you resample,
                          % you MUST scale trialinfo columns 8-15 by
                          % (new_rate/500) before using them as lock samples.
cfg.hp_hz       = 1.0;
cfg.lp_hz       = 100;   % widened vs ds004117's 45 Hz -- gamma needs headroom.
                          % NOTE: the authors' own derivative used 0.5-36 Hz,
                          % which is why we do NOT use data_ica.mat's EEG data,
                          % only its trialinfo.
cfg.line_hz     = 50;    % CONFIRMED via PSD check (50 Hz peak >> 60 Hz peak)

% Channels to drop before preprocessing -- CONFIRMED from the authors'
% Preprocessing_MatlabScript_Fieldtrip.m (their exact exclusion line).
cfg.exclude_channels = [{'M1','M2','EOG'}, ...
    arrayfun(@(i) sprintf('BIP%d',i), 1:24, 'UniformOutput', false)];
% After exclusion, 61 real scalp channels remain (59 listed in channels.tsv
% + Fp1 + Fpz, which channels.tsv omitted -- confirmed by comparing
% channels.tsv against the actual 88-channel raw montage).
cfg.amp_thresh  = 150;   % uV; TODO_CONFIRM once post-preprocessing amplitude
                          % ranges are seen for this dataset specifically.

% ICLabel rejection thresholds for non-Muscle artifact classes (same
% philosophy as Pilot 1's blanket pop_icflag call).
cfg.iclabel.reject_prob = [0.8 1];   % applied to Eye, Heart, LineNoise, ChannelNoise

% Gamma-specific handling for ICLabel "Muscle" components (per the Pilot 1
% decision protocol): a component flagged Muscle may actually be a single
% noisy channel rather than true myogenic contamination. For each
% Muscle-flagged IC (prob in cfg.iclabel.reject_prob range), compute the
% dominant-channel weight from the component topography (EEG.icawinv):
%   z    = zscore(abs(topo))                 -- how extreme the top channel is
%   ratio = top1_weight / top2_weight         -- how dominant vs. the runner-up
% If BOTH exceed threshold -> treat as a single bad channel -> mark that
% channel for interpolation, keep the component (don't remove).
% Otherwise (diffuse multichannel pattern) -> true EMG -> remove the component.
cfg.muscle_ic.prob_thresh   = 0.8;   % ICLabel Muscle probability to consider at all
cfg.muscle_ic.zthresh       = 3.0;   % TODO_CONFIRM: tune once real topographies are seen
cfg.muscle_ic.ratio_thresh  = 3.0;   % TODO_CONFIRM: tune once real topographies are seen

% ---------------- EVENTS ----------------
% Event codes (from ReadMe-EventCodes.txt):
%   Avatar/social cue:    s3021 cue, s3022 objects, s3023 maintenance onset,
%                          s3024 probe, s3025 resp1, s3026 Q2, s3027 resp2
%   Stick/non-social cue: s3041 cue, s3042 objects, s3043 maintenance onset,
%                          s3044 probe, s3045 resp1, s3046 Q2, s3047 resp2
% We do NOT parse these marker pairs ourselves -- see trialinfo note above.
% Kept here only for reference / cross-checking against trialinfo samples.
cfg.ev_maint_avatar = 's3023';
cfg.ev_maint_stick  = 's3043';

% Maintenance window: CONFIRMED empirically (both from raw latency-gap check
% AND independently from trialinfo Maint->LocProbe deltas across trials) to
% be ~1000-1014 ms, essentially fixed with a few ms of jitter. Window chosen
% with margin from both onset (object removal) and offset (probe onset)
% transients.
cfg.maint_window = [0.1 0.9];   % s, relative to trialinfo column 11 (Maint)

% ---------------- trialinfo column map (TrialInformationEEGProcessedData.txt) ----
% 1 MainCon (1-4) | 2 MConCueLR (1-8) | 3 con (1-8) | 4 val (1=valid,2=invalid)
% 5 loc (1=same,2=diff) | 6 cue (1=stick,2=avatar) | 7 LorR (1-4)
% 8 Start | 9 Cue | 10 Targets | 11 Maint | 12 LocProbe | 13 LocRespSamp
% 14 StatQSamp | 15 StatRespSamp | 16 acc_location | 17 acc_status
cfg.ti.col_val        = 4;
cfg.ti.col_loc        = 5;
cfg.ti.col_cue        = 6;   % 1 = stick/nonsocial, 2 = avatar/social
cfg.ti.col_maint      = 11;  % lock sample for maintenance-onset connectivity window
cfg.ti.col_acc_loc    = 16;
cfg.ti.col_acc_stat   = 17;
cfg.ti.cue_stick      = 1;
cfg.ti.cue_avatar     = 2;
cfg.ti.val_valid      = 1;
cfg.ti.val_invalid    = 2;

% Response-time sample columns (for RT-based trial inclusion, see below).
cfg.ti.col_locprobe   = 12;  % location probe onset
cfg.ti.col_locresp    = 13;  % location response
cfg.ti.col_statq      = 14;  % status question onset
cfg.ti.col_statresp   = 15;  % status response

% Trial inclusion: CONFIRMED cfg.ti.col_acc_loc / col_acc_stat are binary
% {0,1} (0 = incorrect, 1 = correct) -- verified via unique() on a real
% trialinfo matrix.
cfg.correct_only     = true;
cfg.accuracy_column  = cfg.ti.col_acc_loc;   % primary: location-probe accuracy
                                              % (most directly relevant to a
                                              % frontoparietal spatial-attention
                                              % analysis); status accuracy is
                                              % available via col_acc_stat if
                                              % you want a secondary check.
cfg.accuracy_correct_value = 1;              % CONFIRMED: 1 = correct, 0 = incorrect

% RT-based inclusion, IN ADDITION to accuracy (per your request: use both).
% RT = response sample - probe/question sample, converted to seconds using
% cfg.resample_hz. Two independent RT filters, one per response:
%   RT_loc  = (col_locresp  - col_locprobe) / srate
%   RT_stat = (col_statresp - col_statq)    / srate
% A trial is excluded if EITHER RT is below the anticipatory floor or is a
% per-subject outlier. Thresholds below are a reasonable starting default,
% NOT a verified fact about this dataset -- adjust after looking at the
% actual RT distributions (s02 should print them).
cfg.rt.use_rt_filter = true;
cfg.rt.min_rt_s      = 0.2;   % TODO_CONFIRM: trials faster than this are
                               % treated as anticipatory/non-genuine responses
cfg.rt.outlier_method = 'mad';  % per-subject: exclude |RT - median| > n*MAD
cfg.rt.outlier_n       = 3;     % TODO_CONFIRM: MAD multiplier

% Degenerate-RT detection: the task paper (Gregory, Wang & Kessler 2022)
% states explicitly there was "no response-window cut-off" -- responses
% were self-paced. So if a subject's RT distribution for either response
% type has almost no spread (MAD below this threshold), that is NOT
% expected task behaviour -- it's a recording/logging anomaly (as found
% for sub-05's RT_stat, confirmed against raw trialinfo samples). In that
% case the RT filter for that specific response type is disabled for that
% subject only (accuracy + the OTHER RT type still apply), and it is
% logged to rt_qc_flags.csv rather than silently trusted or silently
% dropped.
cfg.rt.degenerate_mad_thresh_s = 0.01;   % MAD below 10 ms -> flag as degenerate

cfg.log_power = true;

% ---------------- trialinfo: absolute onset-sample columns (cols 8-10) ---
cfg.ti.col_start          = 8;   % trial start = fixation/ITI onset, ~1s BEFORE cue (NOT t=0)
cfg.ti.col_cue_sample     = 9;   % Cue onset sample -- CONFIRMED t=0 reference (paper's Methods, verbatim)
cfg.ti.col_targets_sample = 10;  % Targets onset sample (~3.5s post-cue per paper Fig. 1 / our data)

% ---------------- ALPHA VALIDATION (s08/s09) ----------------
cfg.alpha.lock_event   = 'cue';         % REVERTED to 'cue' (v2 was right; v3's 'start' was a mistake
                                          % -- see CHANGELOG_alpha_timing_v3v4.md for the full investigation)
cfg.alpha.baseline_s   = [-0.5 -0.1];   % pre-cue baseline, matches paper's Methods EXACTLY ("-500ms to -100ms pre cue
                                          % onset") -- a single fixed window is correct now, no per-trial computation needed.
cfg.alpha.analysis_s   = [-0.5 6.5];    % epoch window rel. to Cue onset: covers baseline through the paper's probe/
                                          % retrieval window (ends at 6s) with a small buffer for wavelet edge effects.
cfg.alpha.foi          = 2:1:30;         % Hz, matches paper
cfg.alpha.cycles       = 3;              % fixed cycles per time window (Morlet), matches paper
cfg.alpha.toi_step_s   = 0.05;           % 50 ms steps, matches paper
cfg.alpha.band         = [8 12];         % alpha band, matches paper
cfg.alpha.avgoverfreq  = true;           % paper: "averaging across the frequency band" BEFORE cluster stats
                                          % (chan x time clustering, not chan x freq x time) -- v3 omitted this.
cfg.alpha.n_perm       = 5000;           % cluster-based permutation count, matches paper exactly (5000 permutations)

% Two a-priori analysis windows, NOW CORRECTLY relative to cue onset
% (paper Fig. 1 / Methods, verbatim):
cfg.alpha.stat_latency_encoding  = [3.5 4.0];  % PRIMARY target: paper found a significant validity effect here (500ms window)
cfg.alpha.stat_latency_retrieval = [5.0 6.0];  % SECONDARY/confirmatory: paper found NO validity effect here (all p>0.1)

% ---------------- FieldTrip path (needed only by s08/s09) ----------------
cfg.fieldtrip_dir = fullfile('..','..','fieldtrip');  % EDIT to your FieldTrip install


% ---------------- ROIs (confirmed present in the real 88-channel montage,
% after excluding M1/M2/EOG/BIP1-24) --------------------------------------
% Per-subject manual bad-channel overrides -- for channels that automated
% clean_rawdata detection missed (e.g. artifact confined to only part of
% the recording, so a whole-recording channel-quality criterion doesn't
% flag it). Found via targeted std() inspection after an s03 sanity check
% flagged abnormal ROI power for a specific subject/condition.
cfg.manual_bad_channels = containers.Map();
cfg.manual_bad_channels('sub-48') = {'P1'};   % std ~153 uV vs ~5-7 uV neighbors;
                                                % power inflation confined to the
                                                % stick-cue block only, confirming
                                                % a block-localized channel artifact
                                                % clean_rawdata's whole-recording
                                                % criterion did not catch.

cfg.roi_frontal  = {'F3','F1','Fz','F2','F4','FC1','FCz','FC2'};
cfg.roi_parietal = {'P3','P1','Pz','P2','P4','CP1','CP2','POz'};  % CPz does NOT

% Path to the FieldTrip neighbour structure (built via ft_prepare_neighbours
% from the montage layout). Consolidated here as the single source of truth --
% s06c, s14, s19, and s20 all previously hardcoded this path independently
% (three of the four to a personal machine path); they now all read
% cfg.neighbours_file instead.
cfg.neighbours_file = fullfile(cfg.deriv_dir, 'neighbours.mat');   % EDIT if not colocated with conn_results.mat
                                                                     % exist on this
                                                                     % montage -- swapped
                                                                     % for POz, confirmed present.

% ---------------- BANDS ----------------
cfg.bands.gamma_low  = [30 45];
cfg.bands.gamma_high = [55 80];

% ---------------- CONNECTIVITY ----------------
cfg.filt_order_factor = 3;

% ---------------- STATS ----------------
cfg.stat.n_perm = 5000;
cfg.stat.alpha  = 0.05;
cfg.stat.tail   = 0;

% ---------------- FIGURES ----------------
cfg.fig.dpi = 300;
end
