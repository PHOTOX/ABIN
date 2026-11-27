#!/usr/bin/env python3
"""Check selected generated inputs, then replay recorded chemistry exports.

Coverage includes BAGEL state/pair requests and orbital reuse, MNDO CI records,
and MOLPRO derivative cards and accuracies. Full input snapshots are checked by
INTERFACES. The supplied geometry is used to check input shape only.
"""

import json
from itertools import combinations
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def main():
    fixture = Path(os.environ["EH_FIXTURE"])
    work = Path.cwd().parent
    launches = work / "chemistry-launches"
    history = launches.read_text().splitlines() if launches.exists() else []
    nstate = int(os.environ["EH_NSTATE"])
    pairs = list(combinations(range(nstate), 2))
    accuracy = int(os.environ["EH_ACCURACY"])
    stage = "gradients"
    if Path.cwd().name == "BAGEL-EH":
        assert sys.argv[1:] == ["input001.json"]
        molecule, reference, derivatives, _ = json.loads(Path(sys.argv[1]).read_text())["bagel"]
        assert molecule["geometry"] and len(molecule["geometry"]) == len(
            (work / "geom.dat.001").read_text().splitlines()
        )
        assert (reference["title"] == "load_ref") == bool(history)
        assert derivatives["method"][0]["nstate"] == nstate
        assert derivatives["method"][0]["thresh"] == 10.0 ** (-accuracy - 1)
        requests = derivatives["grads"]
        assert [job["target"] for job in requests if job["title"] == "force"] == list(range(nstate))
        assert [(job["target"], job["target2"]) for job in requests if job["title"] == "nacme"] == pairs
    else:
        chemistry_input = sys.stdin.read()
        if Path.cwd().name == "MNDO-EH":
            assert not sys.argv[1:]
            if "icross=7 " in chemistry_input:
                stage = "nac"
            assert f"iroot={nstate} " in chemistry_input
            subprocess.run(
                ["awk", "-f", fixture.parents[1] / "validate-mndo-input.awk"],
                input=chemistry_input, text=True, check=True,
            )
        else:
            assert sys.argv[1:] == ["-n", "2", "-s", "--no-xml-output", "-W", str(Path.cwd() / "scratch001")]
            cards = [line for line in chemistry_input.splitlines() if line.startswith("cpmcscf, ")]
            punch = next(line.replace(" ", "") for line in chemistry_input.splitlines() if line.startswith("PUNCH,"))
            if punch == "PUNCH,input001.nacm.pun,new":
                stage = "nac"
            else:
                assert punch == "PUNCH,input001.pun,new"
            states = re.findall(r"^cpmcscf, grad, (\d+)\.1,", chemistry_input, re.M)
            couplings = re.findall(r"^cpmcscf, nacm, (\d+)\.1, (\d+)\.1,", chemistry_input, re.M)
            if stage == "nac":
                assert not states and [(int(i) - 1, int(j) - 1) for i, j in couplings] == pairs
            else:
                assert not couplings and list(map(int, states)) == list(range(1, nstate + 1))
            if stage == "nac" and history[-1:] == ["retry"]:
                accuracy = int(os.environ["EH_FALLBACK_ACCURACY"])
            assert all(card.endswith(f"accu=1d-{accuracy};") for card in cards)

    request = (work / "eh_state.dat").read_text().split()
    count = nstate if stage == "gradients" else len(pairs)
    assert request == [str(nstate)] + ["1"] * count
    with launches.open("a") as log:
        print(stage, file=log)
        if stage == "nac" and "nac" not in history and os.environ["EH_RETRY"] == "1":
            # Fail the first NAC job to exercise ABIN's lower-accuracy retry.
            print("retry", file=log)
            return 1

    step = fixture / "steps" / ("01-gradients" if stage == "gradients" else "02-nac")
    for export in (step / "exports").rglob("*"):
        if export.is_file():
            name = str(export.relative_to(step / "exports"))
            if Path.cwd().name == "MOLPRO-EH":
                # The standalone fixture uses bead 003. Ehrenfest uses bead 001.
                name = name.replace("003", "001")
            output = Path(name)
            output.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(export, output)
    if (step / "stdout").is_file():
        sys.stdout.write((step / "stdout").read_text())
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        (Path.cwd().parent / "mock.failed").touch()
        raise
