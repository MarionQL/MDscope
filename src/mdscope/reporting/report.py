"""Offline evidence export orchestration, with stable PDF and bundle entry points."""
from __future__ import annotations
from pathlib import Path
import tempfile
from datetime import datetime, timezone
from ..config import ensure_dir, config_hash
from .report_inventory import _summary_inventory, _sample_frame_table, _metric_groups, report_table_inventory
from .figure_inventory import figure_inventory
from .pdf import _build_pdf, _maps_figure
from .bundle import _make_bundle

def generate_report(cfg, root, system_repouts=None, force=False):
    """Generate the PDF and optional ZIP after finalization or via CLI report."""
    options = cfg.get("report", {})
    if not force and not options.get("enabled", False):
        return None
    root = Path(root)
    if not (root / "project_summary").is_dir():
        raise FileNotFoundError(
            f"No finalized project_summary under {root}. Run MDscope first."
        )
    if system_repouts is None:
        system_repouts = {
            system: [
                root / system / str(rep.get("label", rep.get("id", f"rep{i+1}")))
                for i, rep in enumerate(spec.get("replicates", []))
            ]
            for system, spec in cfg["systems"].items()
        }
    out = ensure_dir(root / "project_summary" / "report")
    from ..config import atomic_json
    provenance = {
        'report_config_hash': config_hash(cfg), 'sections': options.get('sections',{}),
        'sampling': _summary_inventory(root,cfg,system_repouts).to_dict('records'),
        'statistics_manifest': str(root/'project_summary/statistics/manifest.json'),
        'statistics_methods': cfg.get('statistics', {}),
        'source_analysis_config_status': 'see original run/checkpoint metadata; current rendering config does not prove original trajectory parameters',
        'figure_source_policy': 'PNG content hash plus source path/size/mtime fingerprint; unchanged legacy PNGs are unverified',
    }
    # Render everything in isolation. A failed PDF or bundle export must not
    # replace a previous successful PDF, ZIP, or its provenance sidecars.
    with tempfile.TemporaryDirectory(dir=out, prefix=".mdscope-report-") as directory:
        stage = Path(directory)
        provenance['generated_at'] = datetime.now(timezone.utc).isoformat()
        atomic_json(stage/'report_provenance.json', provenance)
        figure_provenance = out/'figure_provenance.json'
        if figure_provenance.exists():
            from shutil import copyfile
            copyfile(figure_provenance, stage/figure_provenance.name)
        figure_inventory(root,system_repouts,cfg,include_rejected=True).to_csv(stage/'figure_manifest.csv',index=False)
        report_table_inventory(root,cfg).to_csv(stage/'table_manifest.csv',index=False)
        pdf = stage / "MDscope_evidence_report.pdf"
        content = _build_pdf(pdf, cfg, root, system_repouts, options)
        if pdf.stat().st_size < 1000 or pdf.read_bytes()[:5] != b"%PDF-":
            raise RuntimeError("Generated report failed PDF header/size validation")
        atomic_json(stage/'report_content.json', content)
        bundle = _make_bundle(root, pdf, system_repouts, options, cfg) if options.get("bundle_enabled", True) else None
        previous_bundle = out/'MDscope_review_bundle.zip'
        if bundle is None and previous_bundle.exists():
            from ..core.persistence import archive_output
            archive_output(previous_bundle, root/'.archives/report/bundles',
                           config_hash(cfg)+'.zip')
        for path in stage.iterdir():
            path.replace(out/path.name)
    return {"pdf": str(out/pdf.name), "bundle": str(out/bundle.name) if bundle else None}
