import os
import shutil
import argparse
import importlib

import awkward as ak
import dask
import dask_awkward as dak
import egamma_tnp
import mplhep as hep
import numpy as np
import uproot

from matplotlib import pyplot as plt
from coffea.dataset_tools import preprocess
from egamma_tnp import ElectronTagNProbeFromNTuples
from egamma_tnp.utils.histogramming import save_hists

import egamma_tnp.plot
importlib.reload(egamma_tnp.plot)


def pt_low_threshold_plot_setup(**legend_kwargs):
    plt.xlim(10, 400)
    plt.ylim(0, 1.2)
    plt.xlabel(r"Offline electron $P_T$ [GeV]")
    plt.ylabel(r"Efficiency")
    plt.xscale("log")
    plt.xticks([10, 100], [10, 100])
    plt.xticks(
        [20, 30, 40, 50, 60, 70, 80, 90, 200, 300, 400],
        [20, 30, 40, 50, None, None, None, None, 200, 300, 400],
        minor=True,
    )
    plt.legend(**legend_kwargs) if legend_kwargs else plt.legend()


def pt_high_threshold_plot_setup(**legend_kwargs):
    plt.xlim(10, 400)
    plt.ylim(0, 1.2)
    plt.xlabel(r"Offline electron $P_T$ [GeV]")
    plt.ylabel(r"Efficiency")
    plt.legend(**legend_kwargs) if legend_kwargs else plt.legend()


def eta_plot_setup(**legend_kwargs):
    plt.xlim(-2.5, 2.5)
    plt.ylim(0, 1.2)
    plt.xlabel(r"Offline electron $\eta$")
    plt.ylabel(r"Efficiency")
    plt.legend(**legend_kwargs) if legend_kwargs else plt.legend()


def phi_plot_setup(**legend_kwargs):
    plt.xlim(-3.32, 3.32)
    plt.ylim(0, 1.2)
    plt.xlabel(r"Offline electron $\phi$")
    plt.ylabel(r"Efficiency")
    plt.legend(**legend_kwargs) if legend_kwargs else plt.legend()


def plot_eff_only(
    hpass,
    hfail,
    label,
    plottype,
    figure_path,
    legend_kwargs=None,
    cms_kwargs=None,
    eff_kwargs=None,
    efficiency_label="L1T + HLT Efficiency",
):
    legend_kwargs = legend_kwargs or {}
    cms_kwargs = cms_kwargs or {}
    eff_kwargs = eff_kwargs or {}

    hall = hpass + hfail

    num = hpass.values()
    den = hall.values()

    eff = np.divide(
        num,
        den,
        out=np.zeros_like(num, dtype=float),
        where=den > 0,
    )

    err = np.sqrt(
        np.divide(
            eff * (1.0 - eff),
            den,
            out=np.zeros_like(eff, dtype=float),
            where=den > 0,
        )
    )

    edges = hpass.axes[0].edges
    centers = hpass.axes[0].centers
    xerr = np.diff(edges) / 2.0

    plt.figure()

    plt.errorbar(
        centers,
        eff,
        yerr=err,
        xerr=xerr,
        fmt="o",
        label=label,
        **eff_kwargs,
    )

    if plottype == "pt_low_threshold":
        pt_low_threshold_plot_setup(**legend_kwargs)
    elif plottype == "pt_high_threshold":
        pt_high_threshold_plot_setup(**legend_kwargs)
    elif plottype == "eta":
        eta_plot_setup(**legend_kwargs)
    elif plottype == "phi":
        phi_plot_setup(**legend_kwargs)
    else:
        raise ValueError(f"Unknown plottype: {plottype}")

    hep.cms.label(**cms_kwargs)
    plt.ylabel(efficiency_label)

    os.makedirs(os.path.dirname(figure_path), exist_ok=True)
    plt.savefig(figure_path)
    plt.close()


def get_histograms(path):
    with uproot.open(path) as file:
        hpt_barrel_pass = file["pt/barrel/passing"].to_hist()
        hpt_barrel_fail = file["pt/barrel/failing"].to_hist()

        hpt_endcap_loweta_pass = file["pt/endcap_loweta/passing"].to_hist()
        hpt_endcap_loweta_fail = file["pt/endcap_loweta/failing"].to_hist()

        hpt_endcap_higheta_pass = file["pt/endcap_higheta/passing"].to_hist()
        hpt_endcap_higheta_fail = file["pt/endcap_higheta/failing"].to_hist()

        hpt_combined_pass = (
            hpt_barrel_pass
            + hpt_endcap_loweta_pass
            + hpt_endcap_higheta_pass
        )
        hpt_combined_fail = (
            hpt_barrel_fail
            + hpt_endcap_loweta_fail
            + hpt_endcap_higheta_fail
        )

        heta_entire_pass = file["eta/entire/passing"].to_hist()
        heta_entire_fail = file["eta/entire/failing"].to_hist()

        hphi_entire_pass = file["phi/entire/passing"].to_hist()
        hphi_entire_fail = file["phi/entire/failing"].to_hist()

    return (
        hpt_barrel_pass,
        hpt_barrel_fail,
        hpt_endcap_loweta_pass,
        hpt_endcap_loweta_fail,
        hpt_endcap_higheta_pass,
        hpt_endcap_higheta_fail,
        hpt_combined_pass,
        hpt_combined_fail,
        heta_entire_pass,
        heta_entire_fail,
        hphi_entire_pass,
        hphi_entire_fail,
    )


def runfilter(events):
    dataset = events.metadata["dataset"]
    print(dataset, "no run cut")
    return events


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Standalone single-sample Egamma TnP HLT efficiency plotting script"
    )

    parser.add_argument(
        "--outdir",
        dest="outdir",
        help="Output directory",
        type=str,
        default="./",
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Input ROOT file",
    )

    parser.add_argument(
        "--label",
        type=str,
        default="Sample",
        help="Dataset label used in plots and output paths",
    )

    args = parser.parse_args()

    store_dir = args.outdir
    dataset_name = f"events_{args.label}"

    fileset = {
        dataset_name: {
            "files": {
                args.input: "tnpEleTrig/fitter_tree",
            },
        },
    }

    fileset_available, fileset_updated = preprocess(
        fileset,
        step_size=500_000,
        skip_bad_files=True,
    )

    hlt_paths = {
#        "Ele30": "passHLTEle30WPTightGsfTrackIsoFilter",
#        "Ele32": "passHLTEle32WPTightGsfTrackIsoFilter",
#        "Ele115": "passHLTEle115CaloIdVTGsfTrkIdTGsfDphiFilter",
#        "Ele135": "passHLTEle135CaloIdVTGsfTrkIdTGsfDphiFilter",
#        "Ele23Ele12Leg1": "passHLTEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
#        "Ele23Ele12Leg2": "passHLTEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
#        "DoubleEle33SeededLeg": "passHLTEle33CaloIdLMWPMS2Filter",
#        "DoubleEle33UnseededLeg": "passHLTDiEle33CaloIdLMWPMS2UnseededFilter",
      "L1EG20": "passHLTsingleEG23L1matchFilter",
      "singlePhoton30":  "passHLTEG30L1SingleEG20HEFilterLooseHoverE",
      "Ele23": "passHLTsingleEG23ClusterShapeFilter"
    }

    plateau_cuts = {
        "Ele30": 35,
        "Ele32": 35,
        "Ele115": 120,
        "Ele135": 140,
        "Ele23Ele12Leg1": 25,
        "Ele23Ele12Leg2": 15,
        "DoubleEle33SeededLeg": 35,
        "DoubleEle33UnseededLeg": 35,
        "singlePhoton30": 35,
        "Ele23": 25,
        "L1EG20": 25
    }

    tnp = ElectronTagNProbeFromNTuples(
        fileset_available,
        list(hlt_paths.values()),
        cutbased_id="passingCutBasedTight122XV1",
        extra_filter=runfilter,
    )

    to_compute = {}

    for name, trigger in hlt_paths.items():

        if name in ["Ele115", "Ele135"]:
            egamma_tnp.binning.set(
                "pt_bins",
                [
                    5, 10, 15, 20, 22, 26, 28, 30, 32, 34,
                    36, 38, 40, 45, 50, 60, 80, 100, 105,
                    110, 115, 120, 125, 130, 135, 140, 145,
                    150, 200, 250, 300, 350, 400,
                ],
            )
        else:
            egamma_tnp.binning.set(
                "pt_bins",
                [
                    5, 10, 12, 14, 16, 18, 20, 23, 26, 28,
                    30, 32, 34, 36, 38, 40, 45, 50, 60,
                    80, 100, 150, 250, 400,
                ],
            )

        to_compute[name] = tnp.get_1d_pt_eta_phi_tnp_histograms(
            trigger,
            uproot_options={"allow_read_errors_with_report": True},
            eta_regions_pt={
                "barrel": [0.0, 1.4442],
                "endcap_loweta": [1.566, 2.0],
                "endcap_higheta": [2.0, 2.5],
            },
            plateau_cut=plateau_cuts[name],
        )

    dak.necessary_columns(to_compute)

    out, = dask.compute(to_compute)

    plot_dir = f"{store_dir}/plot/{dataset_name}"
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir, exist_ok=True)

    for name, res in out.items():
        hists, report = res

        for dataset, report_arr in report.items():
            os.makedirs(f"{store_dir}/report/{dataset}", exist_ok=True)
            ak.to_json(
                report_arr,
                f"{store_dir}/report/{dataset}/{name}_report.json",
                num_readability_spaces=1,
                num_indent_spaces=4,
            )

        for dataset, hs in hists.items():
            os.makedirs(f"{store_dir}/report/{dataset}", exist_ok=True)
            save_hists(f"{store_dir}/report/{dataset}/{name}_hists.root", hs)

    hep.style.use("CMS")
    hep.style.use(
        {
            "figure.figsize": (6.4, 4.8),
            "font.size": 14,
            "legend.title_fontsize": 14,
            "savefig.bbox": "tight",
        }
    )
    hep.cms.label("Preliminary", data=True, lumi=100, com=15) # ax can be implicit

    rlabel = "2026 (13.6 TeV)"
    eff_kwargs = {"color": "000000"}

    for trigger_name in hlt_paths.keys():

        threshold = trigger_name

        if threshold in ["Ele32", "Ele30"]:
            suffix = "WPTight_Gsf"
        elif threshold in ["Ele115", "Ele135"]:
            suffix = "CaloIdVT_GsfTrkIdT"
        elif threshold == "DoubleEle33SeededLeg":
            suffix = "CaloIdL_MW\nSeeded leg"
        elif threshold == "DoubleEle33UnseededLeg":
            suffix = "CaloIdL_MW\nUnseeded leg"
        elif threshold == "Ele23Ele12Leg1":
            suffix = "CaloIdL_TrackIdL_IsoVL Leg1"
        elif threshold == "Ele23Ele12Leg2":
            suffix = "CaloIdL_TrackIdL_IsoVL Leg2"
        elif threshold == "singlePhoton30":
            suffix = "HoverELoose_L1SingleEG20"
        elif threshold == "Ele23":
             suffix = "CaloIdVT_GsfTrkIdT"
        elif threshold == "L1EG20":
             suffix = "L1SingleEG20"
        else:
            raise ValueError("Couldn't find proper trigger name")


        plateau_cut = plateau_cuts[threshold]

        filename = threshold
        threshold_for_title = (
            threshold
            .replace("Leg1", "")
            .replace("Leg2", "")
            .replace("SeededLeg", "")
            .replace("UnseededLeg", "")
        )

        plottype = (
            "pt_high_threshold"
            if threshold_for_title in ["Ele115", "Ele135"]
            else "pt_low_threshold"
        )

        title = f"HLT_{threshold_for_title}_{suffix}"

        hist_path = f"{store_dir}/report/{dataset_name}/{filename}_hists.root"

        (
            hpt_barrel_pass,
            hpt_barrel_fail,
            hpt_endcap_loweta_pass,
            hpt_endcap_loweta_fail,
            hpt_endcap_higheta_pass,
            hpt_endcap_higheta_fail,
            hpt_combined_pass,
            hpt_combined_fail,
            heta_entire_pass,
            heta_entire_fail,
            hphi_entire_pass,
            hphi_entire_fail,
        ) = get_histograms(hist_path)

        plot_eff_only(
            hpt_barrel_pass,
            hpt_barrel_fail,
            label=args.label,
            plottype=plottype,
            figure_path=f"{store_dir}/report/{filename}_{args.label}_HLT_eff_barrel_pt.pdf",
            legend_kwargs={"title": f"{title}\n$0.00 < |\eta| < 1.44$"},
            cms_kwargs={"loc": 1, "rlabel": rlabel, "text": "Preliminary", "data": True},
            eff_kwargs=eff_kwargs,
        )

        plot_eff_only(
            hpt_endcap_loweta_pass,
            hpt_endcap_loweta_fail,
            label=args.label,
            plottype=plottype,
            figure_path=f"{store_dir}/report/{filename}_{args.label}_HLT_eff_endcap_loweta_pt.pdf",
            legend_kwargs={"title": f"{title}\n$1.57 < |\eta| < 2.00$"},
            cms_kwargs={"loc": 1, "rlabel": rlabel, "text": "Preliminary", "data": True},
            eff_kwargs=eff_kwargs,
        )

        plot_eff_only(
            hpt_endcap_higheta_pass,
            hpt_endcap_higheta_fail,
            label=args.label,
            plottype=plottype,
            figure_path=f"{store_dir}/report/{filename}_{args.label}_HLT_eff_endcap_higheta_pt.pdf",
            legend_kwargs={"title": f"{title}\n$2.00 < |\eta| < 2.50$"},
            cms_kwargs={"loc": 1, "rlabel": rlabel, "text": "Preliminary", "data": True},
            eff_kwargs=eff_kwargs,
        )

        plot_eff_only(
            hpt_combined_pass,
            hpt_combined_fail,
            label=args.label,
            plottype=plottype,
            figure_path=f"{store_dir}/report/{filename}_{args.label}_HLT_eff_combined_pt.pdf",
            legend_kwargs={
                "title": f"{title}\n$0.00 < |\eta| < 1.44$ or $1.57 < |\eta| < 2.50$"
            },
            cms_kwargs={"loc": 1, "rlabel": rlabel, "text": "Preliminary", "data": True},
            eff_kwargs=eff_kwargs,
        )

        plot_eff_only(
            heta_entire_pass,
            heta_entire_fail,
            label=args.label,
            plottype="eta",
            figure_path=f"{store_dir}/report/{filename}_{args.label}_HLT_eff_eta.pdf",
            legend_kwargs={
                "title": f"{title}\n$0.00 < |\eta| < 2.50$\nProbe electron $P_T> {plateau_cut}$ GeV"
            },
            cms_kwargs={"loc": 1, "rlabel": rlabel, "text": "Preliminary", "data": True},
            eff_kwargs=eff_kwargs,
        )

        plot_eff_only(
            hphi_entire_pass,
            hphi_entire_fail,
            label=args.label,
            plottype="phi",
            figure_path=f"{store_dir}/report/{filename}_{args.label}_HLT_eff_phi.pdf",
            legend_kwargs={
                "title": f"{title}\n$0.00 < |\eta| < 2.50$\nProbe electron $P_T> {plateau_cut}$ GeV"
            },
            cms_kwargs={"loc": 1, "rlabel": rlabel, "text": "Preliminary", "data": True},
            eff_kwargs=eff_kwargs,
        )
