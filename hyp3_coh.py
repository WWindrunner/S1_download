"""Submit / fetch HyP3 burst InSAR jobs for event and baseline coherence.

  python hyp3_coh.py submit <bursts.json> <project>  # skips pairs already submitted
  python hyp3_coh.py fetch  <bursts.json> <project>  # download + unzip finished ones

Per burst: event pair (pre1, co) and baseline pair (pre2, pre1), looks 5x1.
Products go to <project>/{event,baseline}/.
"""
import json
import sys
import zipfile
from pathlib import Path

import hyp3_sdk

LOOKS = "5x1"


def wanted_pairs(bursts):
    for b in bursts:
        if b.get("pre1"):
            yield "event", b["pre1"], b["co"]
        if b.get("pre2"):
            yield "baseline", b["pre2"], b["pre1"]


def submit(bursts, project):
    hyp3 = hyp3_sdk.HyP3()
    done = {tuple(j.job_parameters["granules"]) for j in hyp3.find_jobs(name=project)}
    n = 0
    for kind, ref, sec in wanted_pairs(bursts):
        if (ref, sec) in done:
            continue
        hyp3.submit_insar_isce_burst_job(ref, sec, name=project, looks=LOOKS)
        n += 1
    print(f"submitted {n} jobs; {len(done)} already existed")


def fetch(bursts, project):
    hyp3 = hyp3_sdk.HyP3()
    kind_of = {(r, s): k for k, r, s in wanted_pairs(bursts)}
    pending = failed = ready = 0
    for job in hyp3.find_jobs(name=project):
        kind = kind_of.get(tuple(job.job_parameters["granules"]))
        if kind is None:
            continue
        if not job.complete():
            pending += 1
            continue
        if not job.succeeded():
            failed += 1
            print("FAILED:", job.job_parameters["granules"], job.status_code)
            continue
        out = Path(project) / kind
        out.mkdir(parents=True, exist_ok=True)
        ready += 1
        if all((out / Path(f["filename"]).stem).exists() for f in job.files):
            continue
        for z in job.download_files(out):
            with zipfile.ZipFile(z) as zf:
                zf.extractall(out)
            z.unlink()
    print(f"ready: {ready}  pending: {pending}  failed: {failed}")


if __name__ == "__main__":
    action, bursts_json, project = sys.argv[1:4]
    bursts = json.load(open(bursts_json))
    (submit if action == "submit" else fetch)(bursts, project)
