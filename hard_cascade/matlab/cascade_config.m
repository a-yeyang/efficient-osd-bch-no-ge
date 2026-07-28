function cfg = cascade_config(name)
    % Returns one of the two v3 cascade configs (advisor's 3rd task).
    % Port of hc_src/cascade_v3.py CONFIGS_V3.
    switch name
        case 'cfg1_kp4'
            cfg.name = 'Config1 (KP4)';
            cfg.m_rs = 10; cfg.n_rs = 544; cfg.k_rs = 514;
            cfg.m_bch = 8; cfg.t_bch = 1; cfg.n_bch = 144; cfg.k_bch = 136;
        case 'cfg2_255'
            cfg.name = 'Config2 (255)';
            cfg.m_rs = 8; cfg.n_rs = 255; cfg.k_rs = 239;
            cfg.m_bch = 8; cfg.t_bch = 2; cfg.n_bch = 255; cfg.k_bch = 239;
        otherwise
            error('cascade_config: unknown config name %s', name);
    end
    cfg.t_rs = (cfg.n_rs - cfg.k_rs) / 2;
end
