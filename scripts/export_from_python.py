"""Export the exact freMTPL2 train/test split and claims used by motor-pricing-glm (seed 7) so R fits the same data."""
import sys
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split

REPO = Path.home() / "motor-pricing-glm"
sys.path.insert(0, str(REPO / "src"))
from motorpricing import model as m  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "data"
pol, claims = m.load(REPO / "data")
train, test = train_test_split(pol, test_size=0.2, random_state=7)
pol = pol.assign(split=pol["IDpol"].map({**{i: "train" for i in train["IDpol"]}, **{i: "test" for i in test["IDpol"]}}))
pol.to_csv(OUT / "fremtpl2_policies.csv", index=False)
claims.to_csv(OUT / "fremtpl2_claims.csv", index=False)
print(pol.shape, pol["split"].value_counts().to_dict(), claims.shape)
