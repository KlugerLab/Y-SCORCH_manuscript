"""Read revision/config/local_paths.R from Python, so the Python tools use the same
machine-specific locations as the R code (NAME <- "value" lines only)."""
import os, re
_REPO = os.environ.get("SCORCH_REPO") or os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_F = os.path.join(_REPO, "revision", "config", "local_paths.R")
if not os.path.exists(_F):
    raise SystemExit(f"missing {_F} - copy local_paths.R.example and set the data locations")
P = {m.group(1): m.group(2) for m in re.finditer(r'^\s*([A-Z_]+)\s*<-\s*"([^"]+)"', open(_F).read(), re.M)}
