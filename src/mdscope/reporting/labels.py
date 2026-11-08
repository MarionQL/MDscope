"""Scientific display labels shared by plots and PDF tables, without plotting imports."""
METRIC_LABELS = {
    "rmsd_A": ("RMSD", "A"),
    "radius_of_gyration_A": ("Radius of gyration", "A"),
    "ligand_rmsd_A": ("Ligand RMSD", "A"),
    "ligand_min_distance_A": ("Ligand minimum distance", "A"),
    "ligand_contact_residue_count": ("Ligand-contacting residues", "Count"),
    "membrane_thickness_A": ("Membrane thickness", "A"),
    "area_per_lipid_A2": ("Area per lipid", "A^2"),
    "area_per_lipid_upper_A2": ("Box area / upper-leaflet lipid count", "Å²"),
    "area_per_lipid_lower_A2": ("Box area / lower-leaflet lipid count", "Å²"),
    "membrane_msd_lateral_A2": ("Membrane lateral MSD", "A^2"),
    "total_sasa_nm2": ("Total SASA", "nm^2"),
    "helix_fraction": ("Helix fraction", "Fraction"),
    "sheet_fraction": ("Sheet fraction", "Fraction"),
    "coil_fraction": ("Coil fraction", "Fraction"),
}


def _metric_meta(metric: str):
    if metric in METRIC_LABELS:
        return METRIC_LABELS[metric]
    if metric == "lipid_order_mean":
        return "Mean lipid acyl-chain order", "-S_CH"
    if metric.startswith('sdc4_'):
        import re
        label = metric.removeprefix('sdc4_')
        unit = ''
        for suffix, symbol in (('_residue_pairs','residue pairs'),('_neighbor_residues','residues'),('_lipid_contact_lipids','lipids'),('_nm2','nm²'),('_Ainv','1/Å'),('_A2','Å²'),('_A','Å'),('_deg','degrees'),('_fraction','fraction'),('_count','count')):
            if label.endswith(suffix):
                label, unit = (label, symbol) if suffix in {'_residue_pairs','_neighbor_residues','_lipid_contact_lipids'} else (label[:-len(suffix)], symbol)
                break
        label = label.replace('proline_hydrophobic','VVHPLVPL').replace('acidic_patch','acidic patch')
        label = re.sub(r'PRO([AB])_\d+', r'chain \1', label)
        label = label.replace('steric_occlusion','contact occlusion proxy').replace('_',' ')
        return label[:1].upper()+label[1:], unit
    if metric.startswith("lipid_enrichment_") and metric.endswith("_area_density_ratio"):
        label = metric[len("lipid_enrichment_"):-len("_area_density_ratio")].replace("_", " ")
        return f"{label}: density ratio (protein-excluded planar area model)", "ratio"
    if metric.startswith("lipid_order_") and metric.endswith("_mean"):
        label = metric[len("lipid_order_"):-len("_mean")].replace("_", " ")
        return f"Mean lipid order: {label}", "-S_CH"
    if metric.endswith("_count"):
        return metric.replace("_", " ").title(), "Count"
    if metric.endswith("_A"):
        return metric[:-2].replace("_", " ").title(), "A"
    if metric.endswith("_A2"):
        return metric[:-3].replace("_", " ").title(), "A^2"
    if metric.endswith("_fraction"):
        return metric.replace("_", " ").title(), "Fraction"
    return metric.replace("_", " ").title(), metric
