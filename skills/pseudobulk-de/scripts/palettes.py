"""Project palettes and figure choices the house style leaves open.

Fill each field once per project; every figure script reads this file.
  dotplot defaults:    [fill in: size range, colour scale, dot.min]
  volcano thresholds:  [fill in; one seen pair is Q < 0.1 and |log2FC| > 0.5]
  heatmap anchors:     [fill in; one seen example is -4 / 0 / 4]
  font family:         [fill in: e.g. Helvetica]

Palettes are dicts keyed by biological label, never by position.
The two below are placeholders; the project replaces them with its own labels.
"""

pal_condition = {"healthy": "#4477AA", "disease": "#EE6677"}
pal_cell_type = {"t_cell": "#228833", "b_cell": "#CCBB44", "monocyte": "#AA3377"}


def palette_for(labels, pal):
    """Colours for a sequence of labels; an unknown label raises.

    Failing here means a mislabelled or renamed group can never fall through
    as an unplanned colour.
    """
    missing = sorted({str(label) for label in labels} - pal.keys())
    if missing:
        raise KeyError(
            f"palette_for: labels not in the palette: {missing}; "
            "add them to palettes.py or fix the labels."
        )
    return [pal[str(label)] for label in labels]


if __name__ == "__main__":
    assert palette_for(["disease"], pal_condition) == ["#EE6677"]
    try:
        palette_for(["unknown"], pal_condition)
    except KeyError:
        pass
    else:
        raise SystemExit("palette_for accepted an unknown label")
