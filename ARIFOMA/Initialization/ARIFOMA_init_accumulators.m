function acc = ARIFOMA_init_accumulators()
acc.sniffAnnots = struct('who',{},'t0',{},'t1',{},'f0',{},'f1',{},'fmin',{},'fmax',{});
acc.total_chirps_seen = 0;

acc.all_interference_durations      = cell(1);
acc.all_interference_durations_LOS  = cell(1);
acc.all_interference_durations_REFL = cell(1);
acc.all_interference_dist_LOS       = cell(1);
acc.all_interference_dist_REFL      = cell(1);
acc.all_nIntf_per_chirp_total       = cell(1);
acc.all_nIntf_per_chirp_LOS         = cell(1);
acc.all_nIntf_per_chirp_REFL        = cell(1);
acc.all_nOrthogonalDecisionEpochs   = cell(1);
acc.all_r_hat_sub                   = cell(1);

acc.P_coll_per_n            = 0;
acc.total_interfered_chirps = 0;
acc.max_iter                = 1000;

acc.reaction_log = struct('who',{},'eventType',{},'a',{},'b',{},'t0',{},'t1',{},'dt',{},'t_react',{},'cost',{},'iter',{});
acc.ops_log = struct('t_s',{},'who',{},'action',{},'fmin',{},'fmax',{},'score',{},'reason',{},'note',{});
acc.ops_log_events = struct('who',{},'sid',{},'action',{},'t_s',{},'frame',{},'fmin',{},'fmax',{},'score',{},'reason',{},'note',{});

acc.sniff_history_I_all = struct('t0',{},'t1',{},'fmin',{},'fmax',{},'interference',{},'tag',{});
end
