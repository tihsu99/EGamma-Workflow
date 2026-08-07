import argparse
import logging
from pathlib import Path
from time import perf_counter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)
LOGGER.info("Starting plot_all_filters.py")

LOGGER.info("Importing Awkward and Dask")
import awkward as ak
import dask
import dask_awkward as dak

LOGGER.info("Importing plotting and ROOT I/O packages")
import mplhep as hep
import uproot
import yaml
from coffea.dataset_tools import preprocess
from dask.diagnostics import ProgressBar

LOGGER.info("Importing egamma_tnp")
import egamma_tnp
from egamma_tnp import ElectronTagNProbeFromNTuples
from egamma_tnp.plot import plot_ratio
from egamma_tnp.utils.histogramming import save_hists

LOGGER.info("All Python dependencies imported")


PT_BINS = [
    5, 10, 12, 14, 16, 18, 20, 23, 26, 28, 30, 32, 34, 36, 38, 40, 45,
    50, 60, 80, 100, 150, 250, 400,
]
HIGH_PT_BINS = [
    5, 10, 15, 20, 22, 26, 28, 30, 32, 34, 36, 38, 40, 45, 50, 60, 80,
    100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 200, 250, 300,
    350, 400,
]
ETA_REGIONS = {
    "barrel": [0.0, 1.4442],
    "endcap_loweta": [1.566, 2.0],
    "endcap_higheta": [2.0, 2.5],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare the efficiency of every configured HLT filter."
    )
    parser.add_argument("--reference", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--reference-label", default="Reference")
    parser.add_argument("--target-label", default="Target")
    parser.add_argument("--outdir", default=".")
    parser.add_argument(
        "--config",
        default=str(Path(__file__).with_name("config") / "run3_filter_efficiency.yaml"),
    )
    parser.add_argument("--tree", default="tnpEleTrig/fitter_tree")
    parser.add_argument("--cutbased-id", default="passingCutBasedTight122XV1")
    parser.add_argument("--rlabel", default="2026 (13.6 TeV)")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of local Dask worker threads (default: 1)",
    )
    return parser.parse_args()


def load_hlt_paths(config_path):
    with open(config_path, encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    hlt_paths = config.get("hlt_paths", {})
    if not hlt_paths:
        raise ValueError(f"No hlt_paths found in {config_path}")
    for path_name, settings in hlt_paths.items():
        if not settings.get("filters"):
            raise ValueError(f"No filters configured for {path_name}")
        if "plateau_cut" not in settings:
            raise ValueError(f"No plateau_cut configured for {path_name}")
        unknown_overrides = set(settings.get("plateau_cut_overrides", {})) - set(
            settings["filters"]
        )
        if unknown_overrides:
            raise ValueError(
                f"Plateau overrides reference unknown filters in {path_name}: "
                f"{sorted(unknown_overrides)}"
            )
    return hlt_paths


def no_run_cut(events):
    return events


def get_histograms(path):
    with uproot.open(path) as root_file:
        pt_pass = {
            region: root_file[f"pt/{region}/passing"].to_hist()
            for region in ETA_REGIONS
        }
        pt_fail = {
            region: root_file[f"pt/{region}/failing"].to_hist()
            for region in ETA_REGIONS
        }
        pt_pass["combined"] = (
            pt_pass["barrel"]
            + pt_pass["endcap_loweta"]
            + pt_pass["endcap_higheta"]
        )
        pt_fail["combined"] = (
            pt_fail["barrel"]
            + pt_fail["endcap_loweta"]
            + pt_fail["endcap_higheta"]
        )
        return {
            **{
                f"{region}_pt": (pt_pass[region], pt_fail[region])
                for region in (*ETA_REGIONS, "combined")
            },
            "eta": (
                root_file["eta/entire/passing"].to_hist(),
                root_file["eta/entire/failing"].to_hist(),
            ),
            "phi": (
                root_file["phi/entire/passing"].to_hist(),
                root_file["phi/entire/failing"].to_hist(),
            ),
        }


def get_plateau_cut(hlt_settings, filter_name):
    return hlt_settings.get("plateau_cut_overrides", {}).get(
        filter_name, hlt_settings["plateau_cut"]
    )


def plot_settings(plot_name, hlt_settings, filter_name):
    plateau_cut = get_plateau_cut(hlt_settings, filter_name)
    titles = {
        "barrel_pt": r"$0.00 < |\eta| < 1.44$",
        "endcap_loweta_pt": r"$1.57 < |\eta| < 2.00$",
        "endcap_higheta_pt": r"$2.00 < |\eta| < 2.50$",
        "combined_pt": (
            r"$0.00 < |\eta| < 1.44$ or $1.57 < |\eta| < 2.50$"
        ),
        "eta": r"$0.00 < |\eta| < 2.50$"
        + f"\nProbe electron $P_T > {plateau_cut}$ GeV",
        "phi": r"$0.00 < |\eta| < 2.50$"
        + f"\nProbe electron $P_T > {plateau_cut}$ GeV",
    }
    if plot_name.endswith("_pt"):
        plot_type = "pt_high_threshold" if hlt_settings.get("high_threshold") else "pt_low_threshold"
    else:
        plot_type = plot_name
    return plot_type, titles[plot_name]


def write_plot_pair(
    reference_hists,
    target_hists,
    plot_name,
    output_stem,
    hlt_name,
    filter_name,
    hlt_settings,
    args,
):
    plot_type, selection = plot_settings(plot_name, hlt_settings, filter_name)
    for extension in ("pdf", "png"):
        plot_ratio(
            *reference_hists,
            *target_hists,
            label1=args.reference_label,
            label2=args.target_label,
            denominator_type="failing",
            plottype=plot_type,
            figure_path=f"{output_stem}.{extension}",
            legend_kwargs={"title": f"{hlt_name}\n{filter_name}\n{selection}"},
            cms_kwargs={"loc": 1, "rlabel": args.rlabel},
            eff2_kwargs={"color": "#e42536"},
            efficiency_label="L1T + HLT Efficiency",
            ratio_ymin=0.9,
            ratio_ymax=1.1,
        )


def main():
    args = parse_args()
    if args.reference_label == args.target_label:
        raise ValueError("Reference and target labels must be different")
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")

    LOGGER.info("Loading configuration: %s", args.config)
    hlt_paths = load_hlt_paths(args.config)
    filter_count = sum(len(settings["filters"]) for settings in hlt_paths.values())
    LOGGER.info("Configured %d HLT paths and %d filters", len(hlt_paths), filter_count)

    datasets = {
        args.reference_label: {"files": {args.reference: args.tree}},
        args.target_label: {"files": {args.target: args.tree}},
    }
    LOGGER.info("Preprocessing reference and target ROOT files")
    started = perf_counter()
    fileset, _ = preprocess(datasets, step_size=500_000, skip_bad_files=True)
    LOGGER.info("Preprocessing finished in %.1f seconds", perf_counter() - started)
    all_filters = [
        filter_name
        for settings in hlt_paths.values()
        for filter_name in settings["filters"]
    ]
    tnp = ElectronTagNProbeFromNTuples(
        fileset,
        all_filters,
        cutbased_id=args.cutbased_id,
        extra_filter=no_run_cut,
    )

    LOGGER.info("Building the Dask task graph")
    to_compute = {}
    for hlt_name, settings in hlt_paths.items():
        LOGGER.info("  %s: %d filters", hlt_name, len(settings["filters"]))
        egamma_tnp.binning.set(
            "pt_bins", HIGH_PT_BINS if settings.get("high_threshold") else PT_BINS
        )
        to_compute[hlt_name] = {
            filter_name: tnp.get_1d_pt_eta_phi_tnp_histograms(
                filter_name,
                uproot_options={"allow_read_errors_with_report": True},
                eta_regions_pt=ETA_REGIONS,
                plateau_cut=get_plateau_cut(settings, filter_name),
            )
            for filter_name in settings["filters"]
        }

    dak.necessary_columns(to_compute)
    LOGGER.info("Computing all efficiency histograms with %d worker(s)", args.workers)
    started = perf_counter()
    with ProgressBar():
        (results,) = dask.compute(
            to_compute,
            scheduler="threads",
            num_workers=args.workers,
        )
    LOGGER.info("Histogram computation finished in %.1f seconds", perf_counter() - started)
    output_root = Path(args.outdir)

    for hlt_name, settings in hlt_paths.items():
        LOGGER.info("Writing plots for %s", hlt_name)
        for index, filter_name in enumerate(settings["filters"], start=1):
            LOGGER.info(
                "  [%d/%d] %s",
                index,
                len(settings["filters"]),
                filter_name,
            )
            hists, reports = results[hlt_name][filter_name]
            hist_paths = {}
            for dataset, dataset_hists in hists.items():
                hist_path = output_root / "hists" / dataset / hlt_name / f"{filter_name}.root"
                hist_path.parent.mkdir(parents=True, exist_ok=True)
                save_hists(str(hist_path), dataset_hists)
                hist_paths[dataset] = hist_path
            for dataset, report in reports.items():
                report_path = output_root / "reports" / dataset / hlt_name / f"{filter_name}.json"
                report_path.parent.mkdir(parents=True, exist_ok=True)
                ak.to_json(
                    report,
                    str(report_path),
                    num_readability_spaces=1,
                    num_indent_spaces=4,
                )

            reference = get_histograms(hist_paths[args.reference_label])
            target = get_histograms(hist_paths[args.target_label])
            plot_dir = output_root / "plot" / hlt_name / filter_name
            plot_dir.mkdir(parents=True, exist_ok=True)
            for plot_name in reference:
                write_plot_pair(
                    reference[plot_name],
                    target[plot_name],
                    plot_name,
                    plot_dir / f"efficiency_{plot_name}",
                    hlt_name,
                    filter_name,
                    settings,
                    args,
                )
    LOGGER.info("Finished. Output written under %s", output_root.resolve())


if __name__ == "__main__":
    hep.style.use("CMS")
    hep.style.use(
        {
            "figure.figsize": (6.4, 4.8),
            "font.size": 14,
            "legend.title_fontsize": 10,
            "savefig.bbox": "tight",
        }
    )
    main()
