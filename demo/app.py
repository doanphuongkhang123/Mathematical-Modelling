"""#!/usr/bin/env python3
Run:
    python3 -m demo.app
    # then open http://127.0.0.1:8000 in a browser

or:
    python3 demo/app.py --port 8000
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

# --- Make the project's own modules importable ---------------------------------
DEMO_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(DEMO_DIR)
EQUI_DIR = os.path.join(PROJECT_DIR, "equilibria_stability")
for _p in (PROJECT_DIR, EQUI_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from model import PARAMETER_NAMES, params_from_dict          # noqa: E402
from analysis import compute_all_equilibria, p0_stability    # noqa: E402
from stability import jacobian, charpoly, poly_roots         # noqa: E402
from simulator.rhs import rhs                                 # noqa: E402
from simulator.solvers import rk4                             # noqa: E402

INDEX_HTML = os.path.join(DEMO_DIR, "index.html")


def _params(raw: dict) -> "object":
    """Validate + build a Params object from a raw dict of the 9 parameters."""
    clean = {name: float(raw[name]) for name in PARAMETER_NAMES}
    return params_from_dict(clean)


# --- Local stability via the Jacobian (the paper's eigenvalue analysis) --------
# jacobian / charpoly / poly_roots are imported from equilibria_stability/stability.py
# (single shared source of truth, also used by theorem_validation).  This realises
# the Routh-Hurwitz / Hopf-bifurcation study of Sections 2 and 4.

def stability_at(state, params):
    """Eigenvalues + a plain-language classification of local stability."""
    eigs = poly_roots(charpoly(jacobian(state, params)))
    eig_list = [[e.real, e.imag] for e in eigs]
    max_re = max(e.real for e in eigs)
    has_pos = any(e.real > 1e-7 for e in eigs)
    has_neg = any(e.real < -1e-7 for e in eigs)
    spiral = any(abs(e.imag) > 1e-6 for e in eigs)
    # near-Hopf: a complex conjugate pair sitting close to the imaginary axis
    hopf_near = any(abs(e.imag) > 1e-4 and abs(e.real) < 0.05 for e in eigs)

    if max_re < -1e-7:
        kind = "Stable focus (damped oscillation)" if spiral else "Stable node"
    elif max_re > 1e-7:
        if has_neg and has_pos:
            kind = "Saddle (unstable)"
        else:
            kind = "Unstable focus" if spiral else "Unstable node"
    else:
        kind = "Marginal — on a bifurcation"

    return {
        "eigs": eig_list,
        "max_re": max_re,
        "kind": kind,
        "hopf_near": bool(hopf_near and max_re > -0.05),
    }


def regime_label(params) -> dict:
    """Classify the long-term regime from the analytically computed equilibria.

    Mirrors the paper's hierarchy of equilibria P0 (washout) -> P1 (prey only)
    -> P2 (prey + first predator) -> P3 (full coexistence).  The deciding test
    for washout is the P0 stability threshold f1(1) - D1 from p0_stability().
    """
    eq = compute_all_equilibria(params)
    p0_stable, threshold = p0_stability(params)

    if p0_stable:
        label, level = "Washout (P0 stable)", 0
    elif eq["P3"] is not None:
        label, level = "Coexistence (P3)", 3
    elif eq["P2"] is not None:
        label, level = "Prey + 1st predator (P2)", 2
    elif eq["P1"] is not None:
        label, level = "Prey only (P1)", 1
    else:
        label, level = "Washout", 0

    # --- Trophic-level invasion thresholds (transverse eigenvalues) ----------
    # Each is the growth rate of the next species when it is rare at the lower
    # equilibrium: >0 => that level can invade and persist; <0 => it washes out.
    lam_prey = threshold                                   # f1(1) - D1  (P0 -> P1)
    lam_y = None
    if eq["P1"] is not None:
        x1 = eq["P1"][1]
        lam_y = params.f2.value(x1) - params.D2           # f2(x1) - D2 (P1 -> P2)
    lam_z = None
    if eq["P2"] is not None:
        y2 = eq["P2"][2]
        lam_z = params.f3.value(y2) - params.D3           # f3(y2) - D3 (P2 -> P3)

    # --- Dissipativity band (Theorem 2.1): total biomass S+x+y+z attracted to
    #     [1/Dmax, 1/Dmin] with Dmax/Dmin over {1, D1, D2, D3}. ----------------
    Ds = [1.0, params.D1, params.D2, params.D3]
    Dmax, Dmin = max(Ds), min(Ds)

    # --- Local stability of the active (dominant existing) equilibrium -------
    active = next((n for n in ("P3", "P2", "P1", "P0") if eq[n] is not None), "P0")
    stability = stability_at(eq[active], params)

    return {
        "label": label,
        "level": level,
        "p0_stable": bool(p0_stable),
        "threshold": threshold,  # f1(1) - D1; >0 means prey can invade
        "equilibria": {
            k: (list(v) if v is not None else None) for k, v in eq.items()
        },
        "thresholds": {"prey": lam_prey, "y": lam_y, "z": lam_z},
        "biomass_bounds": {"Dmax": Dmax, "Dmin": Dmin,
                           "lower": 1.0 / Dmax, "upper": 1.0 / Dmin},
        "active": active,
        "stability": stability,
    }


def step_trajectory(state, params, span: float, dt: float, samples: int) -> dict:
    """Integrate forward by `span` time units from `state` and subsample.

    Uses the project's RK4 solver and ODE right-hand side directly so the live
    demo shows exactly the dynamics the rest of the project reproduces.
    """
    def rhs_fn(s):
        return rhs(s, params)

    times, states, info = rk4(rhs_fn, state, t_end=span, dt=dt)

    n = len(times)
    if samples >= n:
        idx = range(n)
    else:
        step = (n - 1) / (samples - 1)
        idx = sorted({int(round(i * step)) for i in range(samples)} | {n - 1})

    out_t = [times[i] for i in idx]
    out_s = [states[i] for i in idx]
    return {
        "times": out_t,
        "states": out_s,
        "state": states[-1],          # full-precision endpoint to resume from
        "exploded": bool(info["exploded"]),
    }


# --- HTTP plumbing -------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):
    server_version = "ChemostatLab/1.0"

    def log_message(self, *args):  # keep the console quiet during live polling
        pass

    def _send_json(self, payload: dict, status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_json(self) -> dict:
        length = int(self.headers.get("Content-Length", 0))
        if length <= 0:
            return {}
        return json.loads(self.rfile.read(length).decode("utf-8"))

    def do_GET(self):  # noqa: N802
        if self.path in ("/", "/index.html"):
            try:
                with open(INDEX_HTML, "rb") as f:
                    body = f.read()
            except FileNotFoundError:
                self.send_error(404, "index.html not found")
                return
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        else:
            self.send_error(404, "Not found")

    def do_POST(self):  # noqa: N802
        try:
            data = self._read_json()
            if self.path == "/equilibria":
                params = _params(data["params"])
                self._send_json(regime_label(params))
            elif self.path == "/step":
                params = _params(data["params"])
                state = [float(v) for v in data["state"]]
                span = float(data.get("span", 2.0))
                dt = float(data.get("dt", 0.01))
                samples = int(data.get("samples", 40))
                self._send_json(step_trajectory(state, params, span, dt, samples))
            else:
                self.send_error(404, "Unknown endpoint")
        except (KeyError, ValueError) as exc:
            self._send_json({"error": str(exc)}, status=400)
        except Exception as exc:  # pragma: no cover - defensive
            self._send_json({"error": f"{type(exc).__name__}: {exc}"}, status=500)


def main() -> None:
    parser = argparse.ArgumentParser(description="Food Chain demo server.")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    url = f"http://{args.host}:{args.port}"
    print("=" * 60)
    print(" Interactive Food Chain  (live demo)")
    print("=" * 60)
    print(f" Serving at: {url}")
    print(" Open that address in your browser, then drag the")
    print(" removal-rate knobs (D1, D2, D3) and press Play.")
    print(" Press Ctrl+C to stop.")
    print("=" * 60)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()
