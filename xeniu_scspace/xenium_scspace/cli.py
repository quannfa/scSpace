"""Command-line interface for xeniu-scspace."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="xenium-scspace",
        description="Xenium → scSpace data processing and analysis pipeline",
    )
    sub = parser.add_subparsers(dest="command", required=True, help="Sub-command")

    # ---- prepare ----
    p_prep = sub.add_parser("prepare", help="Prepare scSpace bundle from Xenium data")
    p_prep.add_argument(
        "--from-xenium", type=Path, default=None, metavar="DIR",
        help="Xenium sample output directory (*_outs)",
    )
    p_prep.add_argument("--output-dir", type=Path, default=Path("_l_result/xenium_bundle"), metavar="DIR")
    p_prep.add_argument("--st-type", default="spot", choices=["spot", "image"])
    # Direct CSV mode
    p_prep.add_argument("--sc-data", type=Path, default=None, metavar="FILE")
    p_prep.add_argument("--sc-meta", type=Path, default=None, metavar="FILE")
    p_prep.add_argument("--st-data", type=Path, default=None, metavar="FILE")
    p_prep.add_argument("--st-meta", type=Path, default=None, metavar="FILE")

    # ---- run ----
    p_run = sub.add_parser("run", help="Run scSpace pipeline on a bundle")
    p_run.add_argument("--bundle-dir", type=Path, required=True, metavar="DIR")
    p_run.add_argument("--output-dir", type=Path, default=Path("_l_result/xenium_output"), metavar="DIR")
    p_run.add_argument("--st-type", default="spot", choices=["spot", "image"])
    p_run.add_argument("--n-features", type=int, default=2000)
    p_run.add_argument("--kernel-type", default="primal")
    p_run.add_argument("--dim", type=int, default=50)
    p_run.add_argument("--batch-size", type=int, default=16)
    p_run.add_argument("--hidden-size", type=int, default=128)
    p_run.add_argument("--common-size", type=int, default=2)
    p_run.add_argument("--lr", type=float, default=0.001)
    p_run.add_argument("--epoch-num", type=int, default=1000)
    p_run.add_argument("--Ks", type=int, default=10)
    p_run.add_argument("--Kg", type=int, default=20)
    p_run.add_argument("--res", type=float, default=0.5)
    p_run.add_argument("--target-num", type=int, default=None)
    p_run.add_argument("--random-seed", type=int, default=123)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "prepare":
        if args.from_xenium is not None:
            from .pipeline import prepare_from_xenium

            bundle = prepare_from_xenium(
                xenium_sample_dir=args.from_xenium,
                output_dir=args.output_dir,
                st_type=args.st_type,
            )
            print(f"Bundle written to: {bundle.output_dir}")
            print(f"  sc_data : {bundle.sc_data_path}")
            print(f"  sc_meta : {bundle.sc_meta_path}")
            print(f"  st_data : {bundle.st_data_path}")
            print(f"  st_meta : {bundle.st_meta_path}")
        elif args.sc_data is not None:
            import pandas as pd

            from .pipeline import prepare

            sc_counts = pd.read_csv(args.sc_data, index_col=0)
            sc_meta = pd.read_csv(args.sc_meta, index_col=0)
            st_counts = pd.read_csv(args.st_data, index_col=0)
            st_meta = pd.read_csv(args.st_meta, index_col=0)
            bundle = prepare(sc_counts, sc_meta, st_counts, st_meta, output_dir=args.output_dir)
            print(f"Bundle written to: {bundle.output_dir}")
        else:
            parser.error("Either --from-xenium or --sc-data/--sc-meta/--st-data/--st-meta is required")

    elif args.command == "run":
        from .pipeline import run_scspace_pipeline

        result = run_scspace_pipeline(
            bundle_dir=args.bundle_dir,
            output_dir=args.output_dir,
            st_type=args.st_type,
            n_features=args.n_features,
            kernel_type=args.kernel_type,
            dim=args.dim,
            batch_size=args.batch_size,
            hidden_size=args.hidden_size,
            common_size=args.common_size,
            lr=args.lr,
            epoch_num=args.epoch_num,
            Ks=args.Ks,
            Kg=args.Kg,
            res=args.res,
            target_num=args.target_num,
            random_seed=args.random_seed,
        )
        print(f"Pipeline output: {result['output_dir']}")
        print(f"  pseudo_space : {result['output_dir']}/pseudo_space.csv")
        print(f"  clustering   : {result['output_dir']}/clustering.csv")


if __name__ == "__main__":
    main()