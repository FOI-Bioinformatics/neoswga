"""Do the quantities this tool optimises separate effective panels from the rest?

Plan step 252 of the condition-aware pool design plan. Held out by STUDY, never
by splitting sites within one experiment.

The evidence matrix rules out calibrating coverage or genome recovery: no
dataset here reports breadth at a stated depth. What three of them do carry is a
published per-panel site geometry alongside a binary outcome, so the narrower
question can be asked -- does a panel that scores better on these quantities
tend to be the one that worked?

The predictors are the papers' own published numbers, not quantities recomputed
by this tool. That makes this a check on the metric DEFINITIONS, not on this
implementation of them. It is the stronger reading for a negative result and the
weaker one for a positive.
"""

import itertools
import json
import pathlib

DATA = pathlib.Path("tests/validation/data")
REACH = 3_000

STUDIES = (
    ("clarke_2017_mtb", "most_effective"),
    ("clarke_2017_wolbachia", "effective"),
    ("dwivedi_yu_2023_prevotella", "successful"),
)

# name -> (extractor, higher_is_better)
PREDICTORS = {
    "fg_bg_ratio": (lambda m: m["fg_bg_ratio"], True),
    "fg_gini": (lambda m: m["fg_gini"], False),
    "fg_max_distance": (lambda m: m["fg_max_distance"], False),
    "bg_over_fg_distance": (lambda m: m["bg_mean_distance"] / m["fg_mean_distance"], True),
    "panel_size": (lambda m: m["size"], True),
    # The two baselines the plan names. Reach coverage is a monotone function of
    # mean spacing, so as a RANKING it is the same statistic as site density;
    # both are listed so that equality is visible rather than assumed.
    "site_density": (lambda m: 1.0 / m["fg_mean_distance"], True),
    "reach_coverage": (lambda m: min(1.0, 2.0 * REACH / m["fg_mean_distance"]), True),
}


def load(include_controls):
    out = {}
    for name, label_key in STUDIES:
        data = json.loads((DATA / f"{name}.json").read_text())
        rows = []
        for set_id, entry in data["sets"].items():
            published = entry.get("published")
            if not published:
                continue
            if not include_controls and entry["outcome"].get("negative_control"):
                continue
            rows.append(
                {
                    "id": set_id,
                    "label": bool(entry["outcome"].get(label_key)),
                    "metrics": published,
                }
            )
        out[name] = rows
    return out


def auc(scores, labels, higher_is_better):
    """Rank-based, threshold-free, ties counted as half.

    Written out rather than imported so a reader can see that a tie contributes
    0.5 and not 0 or 1, which matters on samples this small.
    """
    pos = [s for s, y in zip(scores, labels) if y]
    neg = [s for s, y in zip(scores, labels) if not y]
    if not pos or not neg:
        return None
    wins = sum(
        1.0 if (a > b) == higher_is_better and a != b else 0.5 if a == b else 0.0
        for a, b in itertools.product(pos, neg)
    )
    return wins / (len(pos) * len(neg))


def best_threshold(rows, extract, higher_is_better):
    """The cut maximising accuracy on these rows. Ties go to the lower cut."""
    values = sorted({extract(r["metrics"]) for r in rows})
    best, best_acc = None, -1.0
    for cut in values:
        correct = sum(
            ((extract(r["metrics"]) >= cut) == higher_is_better) == r["label"] for r in rows
        )
        if correct / len(rows) > best_acc:
            best, best_acc = cut, correct / len(rows)
    return best, best_acc


def main():
    for include_controls in (True, False):
        tag = "with" if include_controls else "without"
        studies = load(include_controls)
        counts = {k: (len(v), sum(r["label"] for r in v)) for k, v in studies.items()}
        print(f"\n===== {tag} the negative controls =====")
        for name, (n, pos) in counts.items():
            print(f"  {name}: {n} panels, {pos} effective")

        print("\n-- Within-study discrimination (AUC, 0.5 is no information)")
        header = f"{'predictor':<22}" + "".join(f"{n.split('_')[0][:9]:>11}" for n, _ in STUDIES)
        print(header + f"{'mean':>11}")
        for pname, (extract, higher) in PREDICTORS.items():
            values = []
            cells = ""
            for name, _ in STUDIES:
                rows = studies[name]
                a = auc(
                    [extract(r["metrics"]) for r in rows],
                    [r["label"] for r in rows],
                    higher,
                )
                cells += f"{'n/a':>11}" if a is None else f"{a:>11.2f}"
                if a is not None:
                    values.append(a)
            mean = sum(values) / len(values) if values else float("nan")
            print(f"{pname:<22}{cells}{mean:>11.2f}")

        print("\n-- Leave-one-study-out thresholded accuracy")
        print(f"{'predictor':<22}{'held-out acc':>14}{'majority base':>15}")
        for pname, (extract, higher) in PREDICTORS.items():
            correct = total = 0
            for held, _ in STUDIES:
                train = [r for name, _ in STUDIES if name != held for r in studies[name]]
                cut, _ = best_threshold(train, extract, higher)
                for r in studies[held]:
                    predicted = (extract(r["metrics"]) >= cut) == higher
                    correct += predicted == r["label"]
                    total += 1
            everything = [r for name, _ in STUDIES for r in studies[name]]
            majority = max(
                sum(r["label"] for r in everything),
                len(everything) - sum(r["label"] for r in everything),
            ) / len(everything)
            print(f"{pname:<22}{correct / total:>14.2f}{majority:>15.2f}")


main()
