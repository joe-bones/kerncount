"""tumour_volume experiment type (ADR 0002-0004): mCherry tumours in brain sections, one tumour = one n.

    config    [tumour] settings with defaults; declared factors (levels, reference level, display role, source)
    layout    raw experiment folder (one sub-folder per imaging session, flat files) -> layout.csv rows
    segment   one z-stack -> voxel size, threshold, labelled objects, main mass (primary method, ADR 0002)
    measure   per-tumour measurement with a small cache under _work/tumour/
    table     layout + measurements -> the per-tumour readout table (tumours.csv)
    plots     group plots laid out from the factors' display roles
    etype     the registered ExperimentType that chains them

Keep this package importable without numpy / pandas at module level where experiment_types imports it (etype).
"""
