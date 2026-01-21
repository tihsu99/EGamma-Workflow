#!/usr/bin/env python3
import os
import sys
import argparse
import ROOT
from ROOT import TFile, TCanvas, TH1F, TLegend, gStyle, gROOT, TLatex

def setup_root_style():
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadBottomMargin(0.12)

def branch_exists(tree, name):
    return bool(tree.GetBranch(name))

def list_branches(tree):
    out = []
    bl = tree.GetListOfBranches()
    for i in range(bl.GetEntries()):
        out.append(bl.At(i).GetName())
    return out

def should_skip_branch(bname):
    skip_exact = {
        "run", "lumi", "event", "collection_name", "nr_objects"
    }
    if bname in skip_exact:
        return True
    # Often vectors of strings or metadata; skip
    if bname.startswith("nr_") or bname.startswith("collection_"):
        return True
    return False

def default_binning_for(branch_name):
    # Reasonable defaults; you can refine later
    if "eta" in branch_name:
        return 60, -3.0, 3.0
    if "phi" in branch_name:
        return 64, -3.2, 3.2
    if "sigma" in branch_name:
        return 60, -1.0, 10.0
    if "Chi2" in branch_name or "chi2" in branch_name:
        return 60, 0.0, 50.0
    if "Iso" in branch_name or "iso" in branch_name:
        return 60, 0.0, 200.0
    if "et" == branch_name or branch_name.endswith("_et"):
        return 80, 0.0, 3000.0
    if "energy" in branch_name.lower():
        return 80, 0.0, 5000.0
    if "OneOE" in branch_name:
        return 80, -0.5, 0.5
    # Generic fallback
    return 60, 0.0, 1.0

def make_cut(expr, missing=None):
    if missing is None:
        return ""
    return f"({expr} != {missing})"

def draw_hist(tree, expr, hist_name, bins, xmin, xmax, cut=""):
    hist = TH1F(hist_name, "", bins, xmin, xmax)
    hist.Sumw2()
    cmd = f"{expr}>>{hist_name}"
    if cut:
        tree.Draw(cmd, cut, "goff")
    else:
        tree.Draw(cmd, "", "goff")
    return hist

def normalize(hist):
    integ = hist.Integral()
    if integ > 0:
        hist.Scale(1.0 / integ)

def style_hist(hist, color):
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetFillStyle(0)
    hist.SetStats(0)
    hist.GetYaxis().SetTitle("a.u.")

def plot_all(file1_path, file2_path, output_name, legend1, legend2, max_vars=0):
    f1 = TFile(file1_path, "READ")
    f2 = TFile(file2_path, "READ")
    if f1.IsZombie() or f2.IsZombie():
        print("Error: cannot open one of the files")
        return 2

    t1 = f1.Get("egHLTTree")
    t2 = f2.Get("egHLTTree")
    if not t1 or not t2:
        print("Error: egHLTTree not found")
        return 3

    b1 = set(list_branches(t1))
    b2 = set(list_branches(t2))
    common = sorted(list(b1.intersection(b2)))

    # Keep only relevant branches: explicit eg_* and dynamic egvar*
    to_plot = []
    for b in common:
        if should_skip_branch(b):
            continue
        if b.startswith("eg_") or b.startswith("egvar_") or b.startswith("egvarF_") or b.startswith("egvarI_"):
            to_plot.append(b)

    if max_vars > 0:
        to_plot = to_plot[:max_vars]

    print(f"Common branches: {len(common)}")
    print(f"Plotting branches: {len(to_plot)}")

    # Sentinel handling: your dynamic branches usually use -999
    # For explicit branches, you probably did not use sentinel (leave None)
    def sentinel_for(branch_name):
        if branch_name.startswith("egvar_") or branch_name.startswith("egvarF_") or branch_name.startswith("egvarI_"):
            return -999.0
        return None

    for idx, br in enumerate(to_plot):
        bins, xmin, xmax = default_binning_for(br)
        miss = sentinel_for(br)
        cut = make_cut(br, miss)

        h1 = draw_hist(t1, br, f"h1_{idx}", bins, xmin, xmax, cut=cut)
        h2 = draw_hist(t2, br, f"h2_{idx}", bins, xmin, xmax, cut=cut)

        normalize(h1)
        normalize(h2)

        style_hist(h1, ROOT.kBlue)
        style_hist(h2, ROOT.kRed)

        c = TCanvas(f"c_{idx}", br, 800, 1000)

        main_pad = ROOT.TPad(f"main_{idx}", "main", 0.0, 0.3, 1.0, 0.97)
        ratio_pad = ROOT.TPad(f"ratio_{idx}", "ratio", 0.0, 0.0, 1.0, 0.3)
        main_pad.SetLeftMargin(0.12)
        main_pad.SetRightMargin(0.05)
        main_pad.SetTopMargin(0.05)
        main_pad.SetBottomMargin(0.02)
        ratio_pad.SetLeftMargin(0.12)
        ratio_pad.SetRightMargin(0.05)
        ratio_pad.SetTopMargin(0.02)
        ratio_pad.SetBottomMargin(0.25)
        main_pad.Draw()
        ratio_pad.Draw()

        main_pad.cd()
        main_pad.SetLogy(True)
        main_pad.SetGridx(True)
        main_pad.SetGridy(True)

        max_val = max(h1.GetMaximum(), h2.GetMaximum())
        if max_val <= 0:
            max_val = 1.0
        h1.SetMaximum(max_val * 1.5)
        h1.SetMinimum(max_val * 1e-6)

        h1.Draw("hist")
        h2.Draw("hist same")

        leg = TLegend(0.62, 0.75, 0.93, 0.89)
        leg.SetBorderSize(1)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.040)
        leg.AddEntry(h1, legend1, "l")
        leg.AddEntry(h2, legend2, "l")
        leg.Draw()

        tex = TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.035)
        tex.DrawLatexNDC(0.12, 0.97, br)

        ratio_pad.cd()
        ratio_pad.SetGridx(True)
        ratio_pad.SetGridy(True)

        r = h2.Clone(f"r_{idx}")
        r.SetStats(0)
        r.SetLineColor(ROOT.kBlack)
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0.7)
        r.Divide(h1)

        r.SetMinimum(0.0)
        r.SetMaximum(2.0)
        r.GetXaxis().SetTitle(br)
        r.GetXaxis().SetTitleSize(0.10)
        r.GetXaxis().SetLabelSize(0.09)
        r.GetXaxis().SetTitleOffset(0.9)
        r.GetYaxis().SetTitle("Ratio")
        r.GetYaxis().SetTitleSize(0.11)
        r.GetYaxis().SetTitleOffset(0.5)
        r.GetYaxis().SetLabelSize(0.08)

        r.Draw("E1")

        line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()

        safe_name = br.replace(":", "_").replace("/", "_")
        out_png = f"{output_name}_{safe_name}.png"
        c.SaveAs(out_png)
        c.Close()

        h1.Delete()
        h2.Delete()
        r.Delete()

        if (idx + 1) % 20 == 0:
            print(f"Plotted {idx+1}/{len(to_plot)}")

    f1.Close()
    f2.Close()
    print("Done")
    return 0

def main():
    parser = argparse.ArgumentParser(description="Compare all EGamma ntuple branches (explicit + dynamic egvar*).")
    parser.add_argument("file1", help="Reference ROOT file (blue)")
    parser.add_argument("file2", help="Target ROOT file (red)")
    parser.add_argument("output_name", nargs="?", default="egm_all", help="Output base name")
    parser.add_argument("legend1", nargs="?", default="Reference", help="Legend for reference")
    parser.add_argument("legend2", nargs="?", default="Target", help="Legend for target")
    parser.add_argument("--max-vars", type=int, default=0, help="Plot only first N variables (0 = all)")
    args = parser.parse_args()

    if not os.path.exists(args.file1):
        print(f"File not found: {args.file1}")
        return 1
    if not os.path.exists(args.file2):
        print(f"File not found: {args.file2}")
        return 1

    setup_root_style()
    return plot_all(args.file1, args.file2, args.output_name, args.legend1, args.legend2, max_vars=args.max_vars)

if __name__ == "__main__":
    sys.exit(main())

