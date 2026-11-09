"""Audit SDC4 potassium coordination from one CHARMM TPR and optional XTC.

Outputs are per input trajectory. Keep conditions and replicate IDs separate.
Distances and topology masses are in Å and Da, respectively.
"""
from argparse import ArgumentParser
from pathlib import Path

import yaml

from mdscope.mechanisms.potassium_bridge_audit import run_audit


def run(args):
    config = yaml.safe_load(args.config.read_text())
    run_audit(config, args.tpr, args.xtc, args.output, start=args.start,
              stop=args.stop, stride=args.stride, initial_frame=args.initial_frame)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--tpr", required=True, type=Path)
    coordinates = parser.add_mutually_exclusive_group()
    coordinates.add_argument("--xtc", type=Path)
    coordinates.add_argument("--initial-frame", action="store_true",
                             help="Audit the single starting structure stored in the TPR")
    parser.add_argument("--config", type=Path, default=Path("synd4_full_mechanisms.yaml"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--stride", type=int, default=10)
    args = parser.parse_args()
    if args.start < 0 or args.stride < 1 or (args.stop is not None and args.stop <= args.start):
        parser.error("Expected a forward nonnegative frame window and positive stride")
    run(args)
