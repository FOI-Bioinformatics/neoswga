"""Invariants every `per_strain` block must keep, checked in every test that builds one.

Shared rather than asserted once, because the defect it guards was a counting
convention that held for most panels and failed for one shape (a palindromic
oligo, stored under both strand keys at one offset): a single test with an
ordinary panel passed straight over it.
"""


def assert_site_counts_add_up(block):
    """For every measured strain: intact + affected + not assessed == reference sites.

    `reference_sites` is the panel's site count deduplicated over the two
    strands, so this also pins the denominator of every intact fraction. A
    strain that carries no variant must, by the same token, keep every
    assessable site.
    """
    total = block["reference_sites"]
    for name, record in block["per_strain"].items():
        if record["status"] != "measured":
            assert record["intact_sites"]["value"] is None, name
            continue
        intact = record["intact_sites"]["value"]
        affected = record["affected_sites"]["value"]
        not_assessed = record["not_assessed_sites"]["value"]
        assert intact + affected + not_assessed == total, (
            f"{name}: {intact} intact + {affected} affected + {not_assessed} not "
            f"assessed != {total} reference sites"
        )
        proximal = record["affected_three_prime_proximal"]["value"]
        distal = record["affected_distal"]["value"]
        assert proximal + distal == affected, name
        fraction = record["intact_site_fraction"]["value"]
        if fraction is not None:
            assert not_assessed == 0, name
            assert fraction == intact / total, name
        if record["variants"] == 0:
            assert affected == 0, name
            if not_assessed == 0 and total:
                assert fraction == 1.0, name
