"""Small command dispatcher for data audits and experiment planning."""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=[
            "osdr-inventory",
            "osdr-expression",
            "archs4-catalog",
            "split-plan",
            "experiment-plan",
        ],
    )
    args = parser.parse_args(sys.argv[1:2])
    sys.argv = [sys.argv[0], *sys.argv[2:]]
    if args.command == "osdr-inventory":
        from . import osdr_inventory

        osdr_inventory.main()
    elif args.command == "osdr-expression":
        from . import osdr_expression

        osdr_expression.main()
    elif args.command == "archs4-catalog":
        from . import archs4_catalog

        archs4_catalog.main()
    elif args.command == "split-plan":
        from . import split_plan

        split_plan.main()
    else:
        from . import experiment_plan

        experiment_plan.main()


if __name__ == "__main__":
    main()
