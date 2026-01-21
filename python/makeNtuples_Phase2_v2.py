#!/usr/bin/env python3

import argparse
import glob
import os
import subprocess

import ROOT
from DataFormats.FWLite import Events, Handle


ALL_KEYS = [
    "hltEgammaBestGsfTrackVarsUnseeded_Chi2",
    "hltEgammaBestGsfTrackVarsUnseeded_Deta",
    "hltEgammaBestGsfTrackVarsUnseeded_DetaSeed",
    "hltEgammaBestGsfTrackVarsUnseeded_Dphi",
    "hltEgammaBestGsfTrackVarsUnseeded_MissingHits",
    "hltEgammaBestGsfTrackVarsUnseeded_NLayerIT",
    "hltEgammaBestGsfTrackVarsUnseeded_OneOESeedMinusOneOP",
    "hltEgammaBestGsfTrackVarsUnseeded_OneOESuperMinusOneOP",
    "hltEgammaBestGsfTrackVarsUnseeded_ValidHits",
    "hltEgammaBestGsfTrackVarsUnseeded_fbrem",
    "hltEgammaClusterShapeUnseeded",
    "hltEgammaClusterShapeUnseeded_e2x2",
    "hltEgammaClusterShapeUnseeded_sMajor",
    "hltEgammaClusterShapeUnseeded_sMinor",
    "hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5",
    "hltEgammaClusterShapeUnseeded_sigmaIEtaIEta5x5NoiseCleaned",
    "hltEgammaClusterShapeUnseeded_sigmaIPhiIPhi",
    "hltEgammaClusterShapeUnseeded_sigmaIPhiIPhi5x5",
    "hltEgammaClusterShapeUnseeded_sigmaIPhiIPhi5x5NoiseCleaned",
    "hltEgammaEcalPFClusterIsoUnseeded",
    "hltEgammaEleGsfTrackIsoUnseeded",
    "hltEgammaEleL1TrkIsoUnseeded",
    "hltEgammaGsfTrackVarsUnseeded_Chi2",
    "hltEgammaGsfTrackVarsUnseeded_Deta",
    "hltEgammaGsfTrackVarsUnseeded_DetaSeed",
    "hltEgammaGsfTrackVarsUnseeded_Dphi",
    "hltEgammaGsfTrackVarsUnseeded_MissingHits",
    "hltEgammaGsfTrackVarsUnseeded_NLayerIT",
    "hltEgammaGsfTrackVarsUnseeded_OneOESeedMinusOneOP",
    "hltEgammaGsfTrackVarsUnseeded_OneOESuperMinusOneOP",
    "hltEgammaGsfTrackVarsUnseeded_ValidHits",
    "hltEgammaGsfTrackVarsUnseeded_fbrem",
    "hltEgammaHGCALIDVarsUnseeded_hForHOverE",
    "hltEgammaHGCALIDVarsUnseeded_rVar",
    "hltEgammaHGCALIDVarsUnseeded_sigma2uu",
    "hltEgammaHGCALIDVarsUnseeded_sigma2vv",
    "hltEgammaHGCALIDVarsUnseeded_sigma2ww",
    "hltEgammaHGCALIDVarsUnseeded_sigma2xx",
    "hltEgammaHGCALIDVarsUnseeded_sigma2xy",
    "hltEgammaHGCALIDVarsUnseeded_sigma2yy",
    "hltEgammaHGCALIDVarsUnseeded_sigma2yz",
    "hltEgammaHGCALIDVarsUnseeded_sigma2zx",
    "hltEgammaHGCALIDVarsUnseeded_sigma2zz",
    "hltEgammaHGCalLayerClusterIsoUnseeded",
    "hltEgammaHGCalLayerClusterIsoUnseeded_em",
    "hltEgammaHGCalLayerClusterIsoUnseeded_had",
    "hltEgammaHcalPFClusterIsoUnseeded",
    "hltEgammaHollowTrackIsoUnseeded",
    "hltEgammaHoverEUnseeded",
    "hltEgammaPixelMatchVarsUnseeded_s2",
    "hltEgammaR9Unseeded",
    "hltEgammaR9Unseeded_r95x5",
]


def _expand_eos_ls(pattern):
    try:
        res = subprocess.run(
            f"ls {pattern}",
            capture_output=True,
            text=True,
            shell=True,
            check=False,
        )
        if res.returncode == 0 and res.stdout.strip():
            return [x for x in res.stdout.strip().split("\n") if x.strip()]
    except Exception:
        pass
    return []


def resolve_input_files(args):
    files = []
    if args.input_file:
        files = [args.input_file]
    elif args.input_dir:
        pat = os.path.join(args.input_dir, "*.root")
        if "eos" in args.input_dir:
            files = _expand_eos_ls(pat)
        if not files:
            files = glob.glob(pat)
    elif args.input_pattern:
        if "eos" in args.input_pattern:
            files = _expand_eos_ls(args.input_pattern)
        if not files:
            files = glob.glob(args.input_pattern)

    files = [f for f in files if f and f.strip()]
    files.sort()
    return files


def safe_get_seed_info(obj):
    # Returns (rawId, det) or (0, -1) if unavailable
    try:
        sc = obj.superCluster()
        seed = sc.seed().seed()
        return int(seed.rawId()), int(seed.det())
    except Exception:
        return 0, -1


def safe_get_sc_info(obj):
    # Returns (nClusters, rawEnergy, phiWidth) with defaults if unavailable
    try:
        sc = obj.superCluster()
        nclus = int(sc.clusters().size())
        rawE = float(sc.rawEnergy())
        phiW = float(sc.phiWidth())
        return nclus, rawE, phiW
    except Exception:
        return 0, -999.0, -999.0


def main():
    parser = argparse.ArgumentParser(
        description="Create ntuple for Phase2 EGamma trigger objects and store all embedded obj.var keys."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("input_file", nargs="?", help="Input ROOT file")
    input_group.add_argument("--input-dir", "-d", help="Directory with ROOT files")
    input_group.add_argument("--input-pattern", "-p", help="Glob pattern for ROOT files")

    parser.add_argument("--output", "-o", required=True, help="Output ROOT file")
    parser.add_argument("--max-events", "-n", type=int, default=10, help="Max events (default 10)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose")

    args = parser.parse_args()

    input_files = resolve_input_files(args)
    if not input_files:
        raise RuntimeError("No input files found")

    if args.verbose:
        print(f"Input files: {len(input_files)}")
        for f in input_files:
            print(f"  {f}")
        print(f"Output: {args.output}")
        print(f"Max events: {args.max_events}")

    # Product: vector<trigger::EgammaObject> "hltEgammaHLTExtra" "Unseeded" "HLTX"
    eg_handle = Handle("std::vector<trigger::EgammaObject>")
    eg_label = ("hltEgammaHLTExtra", "Unseeded", "HLTX")

    out_file = ROOT.TFile(args.output, "RECREATE")
    tree = ROOT.TTree("egHLTTree", "EGamma HLT Tree")

    # Event info (store once per event, but as vectors for simplicity)
    run = ROOT.std.vector("int")()
    lumi = ROOT.std.vector("int")()
    event = ROOT.std.vector("int")()

    # Per-object basic kinematics and SC seed info
    eg_et = ROOT.std.vector("float")()
    eg_energy = ROOT.std.vector("float")()
    eg_eta = ROOT.std.vector("float")()
    eg_phi = ROOT.std.vector("float")()

    eg_nrClus = ROOT.std.vector("int")()
    eg_rawEnergy = ROOT.std.vector("float")()
    eg_phiWidth = ROOT.std.vector("float")()
    eg_seedId = ROOT.std.vector("unsigned int")()
    eg_seedDet = ROOT.std.vector("int")()

    # Collection info
    collection_name = ROOT.std.vector("string")()
    nr_objects = ROOT.std.vector("int")()

    tree.Branch("run", run)
    tree.Branch("lumi", lumi)
    tree.Branch("event", event)

    tree.Branch("eg_et", eg_et)
    tree.Branch("eg_energy", eg_energy)
    tree.Branch("eg_eta", eg_eta)
    tree.Branch("eg_phi", eg_phi)

    tree.Branch("eg_nrClus", eg_nrClus)
    tree.Branch("eg_rawEnergy", eg_rawEnergy)
    tree.Branch("eg_phiWidth", eg_phiWidth)
    tree.Branch("eg_seedId", eg_seedId)
    tree.Branch("eg_seedDet", eg_seedDet)

    tree.Branch("collection_name", collection_name)
    tree.Branch("nr_objects", nr_objects)

    # Dynamic branches: one vector<float> per var key
    eg_vars = {}
    for k in ALL_KEYS:
        bname = "egvar_" + k
        eg_vars[k] = ROOT.std.vector("float")()
        tree.Branch(bname, eg_vars[k])

    missing = 0

    events = Events(input_files)

    for iev, evt in enumerate(events):
        if iev >= args.max_events:
            break

        # Clear per-event vectors
        run.clear()
        lumi.clear()
        event.clear()

        eg_et.clear()
        eg_energy.clear()
        eg_eta.clear()
        eg_phi.clear()

        eg_nrClus.clear()
        eg_rawEnergy.clear()
        eg_phiWidth.clear()
        eg_seedId.clear()
        eg_seedDet.clear()

        collection_name.clear()
        nr_objects.clear()

        for v in eg_vars.values():
            v.clear()

        # Fill event id
        run.push_back(int(evt.eventAuxiliary().run()))
        lumi.push_back(int(evt.eventAuxiliary().luminosityBlock()))
        event.push_back(int(evt.eventAuxiliary().event()))

        # Read EGamma objects
        evt.getByLabel(eg_label, eg_handle)
        if not eg_handle.isValid():
            if args.verbose:
                print(f"[event {iev}] invalid handle for {eg_label}")
            continue

        egobjs = eg_handle.product()
        nobj = int(egobjs.size())

        collection_name.push_back("hltEgammaHLTExtra:Unseeded:HLTX")
        nr_objects.push_back(nobj)

        if args.verbose:
            print(f"[event {iev}] nEG = {nobj}")

        for obj in egobjs:
            eg_et.push_back(float(obj.et()))
            eg_energy.push_back(float(obj.energy()))
            eg_eta.push_back(float(obj.eta()))
            eg_phi.push_back(float(obj.phi()))

            nclus, rawE, phiW = safe_get_sc_info(obj)
            sid, sdet = safe_get_seed_info(obj)

            eg_nrClus.push_back(int(nclus))
            eg_rawEnergy.push_back(float(rawE))
            eg_phiWidth.push_back(float(phiW))
            eg_seedId.push_back(int(sid))
            eg_seedDet.push_back(int(sdet))

            for k in ALL_KEYS:
                print(k, obj.var(k, missing))
                eg_vars[k].push_back(float(obj.var(k, missing)))

        tree.Fill()

    out_file.Write()
    out_file.Close()

    print(f"Wrote ntuple: {args.output}")


if __name__ == "__main__":
    main()

