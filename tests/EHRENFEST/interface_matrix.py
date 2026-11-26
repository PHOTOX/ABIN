"""Run two-step Ehrenfest NVE trajectories through the real chemistry adapters.

The mock executable replays fixed exports as the geometry moves. This covers
ABIN's interface calls and propagation with those inputs. Coverage stops short
of live electronic structure, long-time conservation, and thermostat dynamics.
"""

from itertools import product
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIXTURES = HERE.parent / "INTERFACES" / "cases"
NSTEPS = 2


def input_text(adapter, settings):
    keywords = "\n".join(f"{key}={value}," for key, value in settings.items())
    return f"""&general
pot='{adapter.lower()}', mdtype='ehrenfest', irest=0, irandom=123,
nstep={NSTEPS}, dt=0.1, nwrite=1, nwritex=1, nwritef=1, xyz_units='angstrom',
/
&thermostat
inose=0, rem_comvel=.false., rem_comrot=.false., scaleveloc=0,
/
&eh
{keywords}
/
"""


def output_table(work, name, columns):
    rows = [
        line.split()
        for line in (work / name).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert len(rows) == NSTEPS and all(len(row) == columns for row in rows), name
    assert all(math.isfinite(float(value)) for row in rows for value in row), name
    return rows


def snapshot(label, values):
    # Keep both columns of the standard runner's side-by-side diff readable.
    return [f"# {label}"] + [
        " ".join(values[start:start + 10])
        for start in range(0, len(values), 10)
    ]


def run_case(binary, fixture, label, settings, expected_error=None, early_error=False):
    adapter, = (directory.name for directory in (fixture / "setup").iterdir() if directory.is_dir())
    geometry = (fixture / "setup" / "geom.dat.003").read_text()
    atoms = [line.split()[0] for line in geometry.splitlines()]
    nstate = settings["nstate"]
    accuracy = settings["nac_accu1"]
    # MNDO and MOLPRO launch separate NAC jobs. BAGEL serves NACs from its gradient cache.
    retry = accuracy > settings["nac_accu2"] and adapter != "BAGEL-EH"
    with TemporaryDirectory(prefix="abin-eh-") as temp:
        work = Path(temp).resolve()
        shutil.copytree(ROOT / "interfaces" / adapter, work / adapter)
        shutil.copytree(fixture / "setup" / adapter, work / adapter, dirs_exist_ok=True)
        executable = work / "mock-program"
        shutil.copy2(HERE / "mock-chemistry.py", executable)
        (work / "SetEnvironment.sh").write_text(
            'export MNDOEXE="$EH_EXECUTABLE" MOLPROEXE="$EH_EXECUTABLE" NCORES=2\n'
        )
        (work / "input.in").write_text(input_text(adapter, settings))
        (work / "mini.xyz").write_text(f"{len(atoms)}\nFixture geometry (angstrom)\n{geometry}")
        velocities = [f"{atom} {0.0003 if i == 0 else -0.0001} 0.0002 0.0"
                      for i, atom in enumerate(atoms)]
        (work / "velocities.in").write_text(
            f"{len(atoms)}\nFixture velocities (atomic units)\n" + "\n".join(velocities) + "\n"
        )
        env = {key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "ENV"}}
        env.update(
            OMP_NUM_THREADS="1", LC_ALL="C", BAGELEXE=str(executable), EH_EXECUTABLE=str(executable),
            FIXTURE_WORK=str(work), EH_FIXTURE=str(fixture), EH_NSTATE=str(nstate), EH_ACCURACY=str(accuracy),
            EH_FALLBACK_ACCURACY=str(settings["nac_accu2"]), EH_RETRY=str(int(retry)),
        )
        result = subprocess.run(
            [binary, "-x", "mini.xyz", "-v", "velocities.in"],
            cwd=work, env=env, capture_output=True, text=True, timeout=30,
        )
        try:
            log = work / "chemistry-launches"
            assert not (work / "mock.failed").exists(), "Mock chemistry rejected the generated input"
            if expected_error:
                assert result.returncode != 0, f"ABIN accepted input expected to fail with {expected_error}"
                assert expected_error in (work / "ERROR").read_text(), f"Expected {expected_error}"
                assert log.exists() != early_error, "Unexpected executable launch state"
                return [f"# {fixture.name} {label}: {expected_error}"]
            assert result.returncode == 0, "ABIN trajectory failed"
            sequence = ["gradients"] if adapter == "BAGEL-EH" else ["gradients", "nac"]
            expected_launches = sequence * (NSTEPS + 1)
            if retry:
                expected_launches[2:2] = ["retry", "nac"]
                assert result.stderr.count("Retrying Ehrenfest NAC calculation") == 1
            assert log.read_text().splitlines() == expected_launches, "Unexpected chemistry launches"
            # Check every PES row against the adapter's independently stored output.
            pes = output_table(work, "PES.dat", nstate + 1)
            gradients = fixture / "steps" / "01-gradients" / "expected" / "engrad.dat.003"
            energies = [float(line) for line in gradients.read_text().splitlines()[:nstate]]
            for row in pes:
                assert all(math.isclose(float(actual), expected, rel_tol=1e-10, abs_tol=1e-10)
                           for actual, expected in zip(row[1:], energies)), "Unexpected state energies"
            # Compare final-frame values with matrix.dat.ref and check table shape at every step.
            population = output_table(work, "pop.dat", nstate + 3)[-1]
            coupling = output_table(work, "dotprod.dat", 1 + nstate * (nstate - 1) // 2)[-1]
            energy = output_table(work, "energies.dat", 4)[-1]
            forces = (work / "forces.xyz").read_text().splitlines()
            assert len(forces) == NSTEPS * (len(atoms) + 2), "Unexpected force trajectory length"
            vectors = [line.split() for line in forces[-len(atoms):]]
            assert [row[0] for row in vectors] == atoms and all(len(row) == 4 for row in vectors)
            components = [value for row in vectors for value in row[1:]]
            assert all(math.isfinite(float(value)) for value in components), "Non-finite force"
            return [
                f"# {fixture.name} {label}",
                *snapshot("dominant state, populations, norm", population[1:]),
                *snapshot("NAC velocity products", coupling[1:]),
                *snapshot("potential, kinetic, total energy", energy[1:]),
                *snapshot("Cartesian forces in atom order", components),
            ]
        except (AssertionError, OSError, ValueError) as error:
            for log in sorted((work / adapter).glob("*out*")):
                if log.is_file():
                    print(f"{log.name}:\n{log.read_text(errors='replace')[-1500:]}", file=sys.stderr)
            raise RuntimeError(f"{fixture.name} {label}: {error}\n{result.stderr}") from error


def main(binary):
    records = []
    for fixture in sorted(FIXTURES.iterdir()):
        # Single-state adapter fixtures have no NAC request and run in INTERFACES.
        nac_request = fixture / "steps" / "02-nac" / "request"
        if not nac_request.is_file():
            continue
        nstate = int(nac_request.read_text().split()[0])
        common = dict(nstate=nstate, couplings="'analytic'")
        # Pair each integrator with two fixed profiles. Starting state, substeps,
        # norm tolerance, and NAC accuracies vary together within these profiles.
        profiles = {
            "ground": dict(istate_init=1, substep=25, popsumthr=0.001, nac_accu1=5, nac_accu2=5),
            "excited": dict(istate_init=nstate, substep=50, popsumthr=0.02, nac_accu1=7, nac_accu2=4),
        }
        for integrator, (profile, settings) in product(("butcher", "rk4", "euler"), profiles.items()):
            records.extend(run_case(binary, fixture, f"{integrator}/{profile}",
                                    dict(common, **settings, integ=f"'{integrator}'")))

        # Vary one invalid &eh value at a time and check rejection before chemistry
        # starts. The strict norm case below instead exercises rejection during propagation.
        # TODO: Add cases for the general-input compatibility checks and iknow override.
        baseline = dict(common, **profiles["ground"], integ="'euler'")
        rejections = (
            ("nstate", 1, "nstate must be >= 2"),
            ("istate_init", nstate + 1, "istate_init is out of bounds"),
            ("substep", 0, "substep must be >= 1"),
            ("integ", "'unknown'", 'Invalid Ehrenfest integrator "unknown"'),
            ("couplings", "'none'", "Ehrenfest requires analytic NAC vectors"),
            ("nac_accu1", 0, "Invalid Ehrenfest NAC accuracy"),
            ("nac_accu2", 0, "Invalid Ehrenfest NAC accuracy"),
        )
        for keyword, value, diagnostic in rejections:
            records.extend(run_case(binary, fixture, f"reject/{keyword}",
                                    dict(baseline, **{keyword: value}), diagnostic, early_error=True))
        strict_norm = dict(common, **profiles["excited"], integ="'euler'")
        strict_norm["popsumthr"] = 1e-14
        records.extend(run_case(binary, fixture, "reject/popsumthr", strict_norm,
                                "Inaccurate Ehrenfest integration led to unphysical populations"))

    if not records:
        raise RuntimeError("No multistate interface fixtures")
    Path("matrix.dat").write_text("\n".join(records) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
