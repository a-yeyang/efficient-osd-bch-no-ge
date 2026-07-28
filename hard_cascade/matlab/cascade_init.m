function codec = cascade_init(cfg, bch_decoder)
    % Build the RS(outer)+BCH(inner) cascade codec for a given config.
    % Port of hc_src/cascade_v3.py CascadeV3Codec.__init__.
    if nargin < 2 || isempty(bch_decoder), bch_decoder = 'direct'; end
    assert(strcmp(bch_decoder, 'conv') || strcmp(bch_decoder, 'direct'), ...
        'cascade_init: bch_decoder must be ''conv'' or ''direct''');
    codec.cfg = cfg;
    codec.bch_decoder = bch_decoder;

    codec.rs = shortened_rs_init(cfg.m_rs, cfg.n_rs, cfg.k_rs);
    codec.bch = bch_init(cfg.m_bch, cfg.t_bch, cfg.n_bch, cfg.k_bch);

    codec.rs_bits = cfg.m_rs * cfg.n_rs;
    if mod(codec.rs_bits, codec.bch.k) ~= 0
        codec.n_pad_bits = codec.bch.k - mod(codec.rs_bits, codec.bch.k);
    else
        codec.n_pad_bits = 0;
    end
    codec.n_bch_blocks = (codec.rs_bits + codec.n_pad_bits) / codec.bch.k;
    codec.n_coded_bits = codec.n_bch_blocks * codec.bch.n;

    codec.n_info_bits = cfg.k_rs * cfg.m_rs;
    codec.effective_rate = codec.n_info_bits / codec.n_coded_bits;
end
