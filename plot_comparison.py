import os
import shutil
import awkward as ak
import dask
import dask_awkward as dak
import mplhep as hep
import uproot
from coffea.dataset_tools import preprocess
from distributed import Client
from matplotlib import pyplot as plt
import argparse
import egamma_tnp
from egamma_tnp import ElectronTagNProbeFromNTuples
from egamma_tnp.plot import plot_ratio
from egamma_tnp.utils.histogramming import save_hists
import importlib, egamma_tnp.plot
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



lumis = {
    "2025G": 0,
}


if __name__ == "__main__":

  usage  = 'usage: %prog [options]'
  parser = argparse.ArgumentParser(description=usage)
  parser.add_argument('--outdir',     dest = 'outdir',     help='output directory',   type=str, default='./')
  parser.add_argument('--reference',  type=str)
  parser.add_argument('--target', type=str)
  parser.add_argument('--reference_label', type=str, default="Reference")
  parser.add_argument('--target_label', type=str, default="Target")
  args = parser.parse_args()

  fileset = {
    f"events_{args.reference_label}" : {
        "files": {
        args.reference: "tnpEleTrig/fitter_tree",
        },
    },
    f"events_{args.target_label}" : {
        "files": {
        args.target: "tnpEleTrig/fitter_tree",
        },
    },
  }
  store_dir = args.outdir
  fileset_available, fileset_updated = preprocess(fileset, step_size=500_000, skip_bad_files=True)

  # Define the HLT paths for wihich the efficiency plots are needed along with the last passing Filter

  hlt_paths = {
    "Ele30": "passHLTEle30WPTightGsfTrackIsoFilter",
    "Ele32": "passHLTEle32WPTightGsfTrackIsoFilter",
    "Ele115": "passHLTEle115CaloIdVTGsfTrkIdTGsfDphiFilter",
    "Ele135": "passHLTEle135CaloIdVTGsfTrkIdTGsfDphiFilter",
    "Ele23Ele12Leg1": "passHLTEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",
    "Ele23Ele12Leg2": "passHLTEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",
    "DoubleEle33SeededLeg": "passHLTEle33CaloIdLMWPMS2Filter",
    "DoubleEle33UnseededLeg": "passHLTDiEle33CaloIdLMWPMS2UnseededFilter",
#    "SinglePhoton30":  "passHLTEG30L1SingleEG20HEFilterLooseHoverE",
#    "Ele23": "passHLTsingleEG23TrackIsoL1matchFilter"
  }

  # Check what the purpose of defining the plateau cuts is

  plateau_cuts = {
    "Ele30": 35,
    "Ele32": 35,
    "Ele115": 120,
    "Ele135": 140,
    "Ele23Ele12Leg1": 25,
    "Ele23Ele12Leg2": 15,
    "DoubleEle33SeededLeg": 35,
    "DoubleEle33UnseededLeg": 35,
    "Ele23": 25,
    "singlePhoton": 35
  }

  def runfilter(events):
    dataset = events.metadata["dataset"]
    print(dataset, "no run cut")
    return events

  tnp = ElectronTagNProbeFromNTuples(fileset_available, list(hlt_paths.values()), cutbased_id="passingCutBasedTight122XV1", extra_filter=runfilter)

  to_compute = {}

  for name, trigger in hlt_paths.items():
    if name == "Ele115" or name == "Ele135":
        egamma_tnp.binning.set(
            "pt_bins",
            [
                5,
                10,
                15,
                20,
                22,
                26,
                28,
                30,
                32,
                34,
                36,
                38,
                40,
                45,
                50,
                60,
                80,
                100,
                105,
                110,
                115,
                120,
                125,
                130,
                135,
                140,
                145,
                150,
                200,
                250,
                300,
                350,
                400,
            ],
        )
    else:
        egamma_tnp.binning.set(
            "pt_bins",
            [
                5,
                10,
                12,
                14,
                16,
                18,
                20,
                23,
                26,
                28,
                30,
                32,
                34,
                36,
                38,
                40,
                45,
                50,
                60,
                80,
                100,
                150,
                250,
                400,
            ],
        )
    plateau_cut = plateau_cuts[name]
    to_compute[name] = tnp.get_1d_pt_eta_phi_tnp_histograms(
        trigger,
        uproot_options={"allow_read_errors_with_report": True},
        eta_regions_pt={
            "barrel": [0.0, 1.4442],
            "endcap_loweta": [1.566, 2.0],
            "endcap_higheta": [2.0, 2.5],
        },
        plateau_cut=plateau_cut,
    )

  dak.necessary_columns(to_compute)

  (out,) = dask.compute(to_compute)

  for dataset in out["Ele30"][1].keys():
    path = f"{store_dir}/plot/{dataset}"
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)

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

  def get_histograms(path):
    with uproot.open(path) as file:
        hpt_barrel_pass = file["pt/barrel/passing"].to_hist()
        hpt_barrel_fail = file["pt/barrel/failing"].to_hist()
        hpt_endcap_loweta_pass = file["pt/endcap_loweta/passing"].to_hist()
        hpt_endcap_loweta_fail = file["pt/endcap_loweta/failing"].to_hist()
        hpt_endcap_higheta_pass = file["pt/endcap_higheta/passing"].to_hist()
        hpt_endcap_higheta_fail = file["pt/endcap_higheta/failing"].to_hist()
        hpt_combined_pass = hpt_barrel_pass + hpt_endcap_loweta_pass + hpt_endcap_higheta_pass
        hpt_combined_fail = hpt_barrel_fail + hpt_endcap_loweta_fail + hpt_endcap_higheta_fail

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
  ratio_min = 0.9
  ratio_max = 1.1
  for tala in list(hlt_paths.keys()):
    for data_period in [f"events_{args.reference_label}"]:
        for mc_dataset in [f"events_{args.target_label}"]: ## data_period = reference, mc_dataset = target
            tocompare = [data_period, mc_dataset]
            run = []
            for folder in tocompare:
                run.append(folder.split("_", 2)[2][3:] if "data" in folder else folder.split("_", 1)[1])
            threshold = tala
            if threshold == "Ele32" or threshold == "Ele30":
                suffix = "WPTight_Gsf"
            elif threshold == "Ele115" or threshold == "Ele135":
                suffix = "CaloIdVT_GsfTrkIdT"
            elif threshold == "DoubleEle33SeededLeg":
                suffix = "CaloIdL_MW\nSeeded leg"
            elif threshold == "DoubleEle33UnseededLeg":
                suffix = "CaloIdL_MW\nUnseeded leg"
            elif threshold == "Ele23Ele12Leg1":
                suffix = "CaloIdL_TrackIdL_IsoVL Leg1"
            elif threshold == "Ele23Ele12Leg2":
                suffix = "CaloIdL_TrackIdL_IsoVL Leg2"
            else:
                raise ValueError("Couldn't find proper trigger name")

            plateau_cut_dict = {
                "Ele30": 35,
                "Ele32": 35,
                "Ele115": 120,
                "Ele135": 140,
                "Ele23Ele12Leg1": 25,
                "Ele23Ele12Leg2": 15,
                "DoubleEle33SeededLeg": 35,
                "DoubleEle33UnseededLeg": 35,
                "singlePhoton": 35,
                "Ele23": 25
            }
            plateau_cut = plateau_cut_dict[threshold]

            filename = threshold
            threshold = threshold.replace("Leg1", "").replace("Leg2", "").replace("SeededLeg", "").replace("UnseededLeg", "")

            plottype = "pt_high_threshold" if threshold == "Ele115" or threshold == "Ele135" else "pt_low_threshold"
            title = f"HLT_{threshold}_{suffix}"

            rlabel = "2026 (13.6 TeV)"

            if mc_dataset == f"events_{args.target_label}":
                eff2_kwargs = {"color": "#e42536"} # Red                      
            else:
                eff2_kwargs = {"color": "#5790fc"}

            (
                hpt_barrel_pass1,
                hpt_barrel_all1,
                hpt_endcap_loweta_pass1,
                hpt_endcap_loweta_all1,
                hpt_endcap_higheta_pass1,
                hpt_endcap_higheta_all1,
                hpt_combined_pass1,
                hpt_combined_all1,
                heta_entire_pass1,
                heta_entire_all1,
                hphi_entire_pass1,
                hphi_entire_all1,
            ) = get_histograms(f"{store_dir}/report/{tocompare[0]}/{filename}_hists.root")

            (
                hpt_barrel_pass2,
                hpt_barrel_all2,
                hpt_endcap_loweta_pass2,
                hpt_endcap_loweta_all2,
                hpt_endcap_higheta_pass2,
                hpt_endcap_higheta_all2,
                hpt_combined_pass2,
                hpt_combined_all2,
                heta_entire_pass2,
                heta_entire_all2,
                hphi_entire_pass2,
                hphi_entire_all2,
            ) = get_histograms(f"{store_dir}/report/{tocompare[1]}/{filename}_hists.root")

            plot_ratio(
                hpt_barrel_pass1,
                hpt_barrel_all1,
                hpt_barrel_pass2,
                hpt_barrel_all2,
                label1=f"{run[0]}",
                label2=f"{run[1]}",
                denominator_type="failing",
                plottype=plottype,
                figure_path=f"{store_dir}/report/{filename}_{run[0]}_vs_{run[1]}_HLT_eff_barrel_pt.pdf",
                legend_kwargs={"title": f"{title}\n$0.00 < |\eta| < 1.44$"},
                cms_kwargs={"loc": 1, "rlabel": rlabel},
                eff2_kwargs=eff2_kwargs,
                efficiency_label="L1T + HLT Efficiency",
                ratio_ymin=ratio_min,
                ratio_ymax=ratio_max
            )

            plot_ratio(
                hpt_endcap_loweta_pass1,
                hpt_endcap_loweta_all1,
                hpt_endcap_loweta_pass2,
                hpt_endcap_loweta_all2,
                label1=f"{run[0]}",
                label2=f"{run[1]}",
                denominator_type="failing",
                plottype=plottype,
                figure_path=f"{store_dir}/report/{filename}_{run[0]}_vs_{run[1]}_HLT_eff_endcap_loweta_pt.pdf",
                legend_kwargs={"title": f"{title}\n$1.57 < |\eta| < 2.00$"},
                cms_kwargs={"loc": 1, "rlabel": rlabel},
                eff2_kwargs=eff2_kwargs,
                efficiency_label="L1T + HLT Efficiency",
                ratio_ymin=ratio_min,
                ratio_ymax=ratio_max
            )

            plot_ratio(
                hpt_endcap_higheta_pass1,
                hpt_endcap_higheta_all1,
                hpt_endcap_higheta_pass2,
                hpt_endcap_higheta_all2,
                label1=f"{run[0]}",
                label2=f"{run[1]}",
                denominator_type="failing",
                plottype=plottype,
                figure_path=f"{store_dir}/report/{filename}_{run[0]}_vs_{run[1]}_HLT_eff_endcap_higheta_pt.pdf",
                legend_kwargs={"title": f"{title}\n$2.00 < |\eta| < 2.50$"},
                cms_kwargs={"loc": 1, "rlabel": rlabel},
                eff2_kwargs=eff2_kwargs,
                efficiency_label="L1T + HLT Efficiency",
                ratio_ymin=ratio_min,
                ratio_ymax=ratio_max
            )

            plot_ratio(
                hpt_combined_pass1,
                hpt_combined_all1,
                hpt_combined_pass2,
                hpt_combined_all2,
                label1=f"{run[0]}",
                label2=f"{run[1]}",
                denominator_type="failing",
                plottype=plottype,
                figure_path=f"{store_dir}/report/{filename}_{run[0]}_vs_{run[1]}_HLT_eff_combined_pt.pdf",
                legend_kwargs={"title": f"{title}\n$0.00 < |\eta| < 1.44$ or $1.57 < |\eta| < 2.50$"},
                cms_kwargs={"loc": 1, "rlabel": rlabel},
                eff2_kwargs=eff2_kwargs,
                efficiency_label="L1T + HLT Efficiency",
                ratio_ymin=ratio_min,
                ratio_ymax=ratio_max
            )

            plot_ratio(
                heta_entire_pass1,
                heta_entire_all1,
                heta_entire_pass2,
                heta_entire_all2,
                label1=f"{run[0]}",
                label2=f"{run[1]}",
                denominator_type="failing",
                plottype="eta",
                figure_path=f"{store_dir}/report/{filename}_{run[0]}_vs_{run[1]}_HLT_eff_eta.pdf",
                legend_kwargs={"title": f"{title}\n$0.00 < |\eta| < 2.50$\nProbe electron $P_T> {plateau_cut}$ GeV"},
                cms_kwargs={"loc": 1, "rlabel": rlabel},
                eff2_kwargs=eff2_kwargs,
                efficiency_label="L1T + HLT Efficiency",
                ratio_ymin=ratio_min,
                ratio_ymax=ratio_max
            )

            plot_ratio(
                hphi_entire_pass1,
                hphi_entire_all1,
                hphi_entire_pass2,
                hphi_entire_all2,
                label1=f"{run[0]}",
                label2=f"{run[1]}",
                denominator_type="failing",
                plottype="phi",
                figure_path=f"{store_dir}/report/{filename}_{run[0]}_vs_{run[1]}_HLT_eff_phi.pdf",
                legend_kwargs={"title": f"{title}\n$0.00 < |\eta| < 2.50$\nProbe electron $P_T> {plateau_cut}$ GeV"},
                cms_kwargs={"loc": 1, "rlabel": rlabel},
                eff2_kwargs=eff2_kwargs,
                efficiency_label="L1T + HLT Efficiency",
                ratio_ymin=ratio_min,
                ratio_ymax=ratio_max
            )
