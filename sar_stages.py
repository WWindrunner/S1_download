"""Per-scene checkpoints: successful outputs are reused only for matching inputs."""
import json
from pathlib import Path
import time

from nisar_common import file_stamp, fingerprint, validate_raster, write_json


class SceneStages:
    def __init__(self, directory, force=False):
        self.path = Path(directory) / "processing_status.json"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            self.state = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            self.state = {"stages": {}}
        self.force = force
        self.state["complete"] = False
        write_json(self.path, self.state)

    def run(self, name, inputs, action, reference=None, require_valid=True):
        key = fingerprint(inputs)
        record = self.state["stages"].get(name, {})
        if not self.force and record.get("key") == key and record.get("status") == "complete":
            try:
                paths = record["outputs"]
                if [file_stamp(p) for p in paths] != record["output_stamps"]:
                    raise ValueError("Output changed")
                for path in paths:
                    validate_raster(path, reference, require_valid)
                print(f"Reusing completed stage: {name}")
                return paths
            except (OSError, ValueError, KeyError):
                pass
        self.state["stages"][name] = dict(status="running", key=key)
        write_json(self.path, self.state)
        started = time.perf_counter()
        try:
            outputs = [str(Path(p).resolve()) for p in action()]
            if not outputs:
                raise ValueError(f"Stage produced no outputs: {name}")
            for path in outputs:
                validate_raster(path, reference, require_valid)
            self.state["stages"][name] = dict(status="complete", key=key, outputs=outputs,
                                              output_stamps=[file_stamp(p) for p in outputs])
        except Exception as exc:
            self.state["stages"][name] = dict(status="failed", key=key, error=type(exc).__name__)
            write_json(self.path, self.state)
            raise
        write_json(self.path, self.state)
        print(f"{name} completed in {time.perf_counter() - started:.1f} seconds")
        return outputs

    def finish(self):
        self.state["complete"] = True
        write_json(self.path, self.state)


def implementation_stamp(*names):
    root = Path(__file__).resolve().parent
    # Include code hashes, so algorithm edits invalidate affected checkpoints.
    import hashlib
    return {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in names}
