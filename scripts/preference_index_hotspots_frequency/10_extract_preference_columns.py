#!/usr/bin/env python3
"""Extract compact preference-index columns from a full PI summary."""
from __future__ import annotations

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--columns", default="Component,Pref Index (Component),Pref Index (Global)")
    args = parser.parse_args()
    cols = [c.strip() for c in args.columns.split(',')]
    df = pd.read_csv(args.input)
    df[cols].to_csv(args.output, index=False)
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
