#!/usr/bin/env sh
# Check that the *current* environment (conda or not) has all required deps,
# regardless of the env's name. If anything is missing or version-mismatched,
# show a summary and ask whether to proceed.

# Usage:
#   ./check_deps.sh "ete3 numpy>=1.24 pandas==2.2.*"
# or set REQ_SPECS before calling.
REQ_SPECS="${1:-${REQ_SPECS:-'ete3 numpy pandas'}}"   ## space/comma-separated; supports ==,!=,>=,>,<=,<

PROMPT_ON_FAIL=1  ## set to 0 to skip the prompt and fail hard on missing deps

## heads-up if no conda env is active
if [ -z "${CONDA_PREFIX:-}" ]; then
  echo "Note: no conda environment active (CONDA_PREFIX unset). Checking current interpreter anyway." >&2
fi

RESULT="$(
  REQ_SPECS="$REQ_SPECS" python - <<'PY'
import os, sys, json, re, subprocess

req = os.environ.get("REQ_SPECS","").strip()

def _list_conda():
    try:
        out = subprocess.check_output(["conda","list","--json"], stderr=subprocess.DEVNULL)
        data = json.loads(out)
        return {d["name"].lower(): d.get("version","") for d in data}
    except Exception:
        return {}

def _list_pip():
    try:
        out = subprocess.check_output([sys.executable,"-m","pip","list","--format=json"], stderr=subprocess.DEVNULL)
        data = json.loads(out)
        return {d["name"].lower(): d.get("version","") for d in data}
    except Exception:
        return {}

installed = {}
installed.update(_list_conda())
installed.update(_list_pip())

specs = []
for tok in re.split(r"[,\s]+", req):
    if not tok: 
        continue
    m = re.match(r"^([A-Za-z0-9_.\-]+)\s*([!<>=]{1,2})?\s*([A-Za-z0-9*_.+\-]+)?$", tok)
    if not m:
        specs.append((tok,"","",False,"invalid spec"))
        continue
    name, op, ver = m.group(1).lower(), m.group(2) or "", m.group(3) or ""
    specs.append((name, op, ver, True, ""))

def cmp_ok(iv, op, rv):
    if not op:
        return True
    try:
        try:
            from packaging.version import Version as V  # preferred
            iv, rv = V(iv), V(rv)
        except Exception:
            from pkg_resources import parse_version as pv  # fallback (setuptools)
            iv, rv = pv(iv), pv(rv)
    except Exception:
        return (op in ("==","=")) and (iv == rv)
    return ((op in ("==","=") and iv == rv) or
            (op == "!=" and iv != rv) or
            (op == ">"  and iv >  rv) or
            (op == ">=" and iv >= rv) or
            (op == "<"  and iv <  rv) or
            (op == "<=" and iv <= rv))

missing = []
for name, op, ver, ok, why in specs:
    if not ok:
        missing.append((name, op, ver, "invalid requirement"))
        continue
    inst = installed.get(name)
    if inst is None:
        missing.append((name, op, ver, "not installed"))
    elif not cmp_ok(inst, op, ver):
        want = f"{name}{op}{ver}" if op or ver else name
        missing.append((name, op, ver, f"installed {inst}"))

if missing:
    print("MISSING")
    for (name, op, ver, why) in missing:
        want = f"{name}{op}{ver}" if op or ver else name
        print(f"{want}\t{why}")
else:
    print("OK")
PY
)"

status=$(printf '%s\n' "$RESULT" | head -n1)
if [ "$status" = "OK" ]; then
  ## all dependencies satisfied
  exit 0
else
  echo "Dependency check failed:" >&2
  printf '%s\n' "$RESULT" | tail -n +2 | sed 's/^/  - /' >&2
  if [ "$PROMPT_ON_FAIL" -eq 1 ] && [ -t 0 ]; then
    printf "Proceed anyway? [y/N] " >&2
    read -r ans
    case "$ans" in
      y|Y|yes|YES) exit 0 ;;
      *) echo "Aborting. Activate/install required packages and retry." >&2; exit 1 ;;
    esac
  else
    exit 1
  fi
fi
