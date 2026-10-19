# ==== Shared functions (cdri_lib): identical copy in every notebook from 01 onwards ====
# Change a function here, rebuild the notebooks, and every notebook uses the same version.
import os, sys, re, json, math, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
pd.set_option('display.width', 160)

KT2MS, NMI2KM, R_EARTH_KM = 0.514444, 1.852, 6371.0
SURF_TO_GRADIENT = 0.85          # surface-to-gradient wind ratio used only to estimate the Holland B parameter
CRS_GEO, CRS_M = 'EPSG:4326', 'EPSG:32645'

# WMO/TD-1555 (Harper, Kepert & Ginger 2010) Table 1.1: 3-s gust factor G(3, To) at +10 m in tropical cyclones
GUST3 = {'in-land':  {60: 1.49, 180: 1.58, 600: 1.66},
         'off-land': {60: 1.36, 180: 1.44, 600: 1.52},
         'off-sea':  {60: 1.23, 180: 1.31, 600: 1.38},
         'at-sea':   {60: 1.11, 180: 1.17, 600: 1.23}}


def setup_project():
    """Mount Google Drive (in Colab) and return the project root and folder dictionary."""
    if 'google.colab' in sys.modules:
        from google.colab import drive
        try:
            drive.mount('/content/drive')
            root = Path('/content/drive/MyDrive/CDRI_Odisha')
        except Exception as e:
            print('Drive mount failed:', e, '\nFix: Runtime → Disconnect and delete runtime, rerun, use the same Google '
                  'account and tick all permission boxes. Using temporary storage for now.')
            root = Path('/content/CDRI_Odisha')
    else:
        root = Path(os.environ.get('CDRI_ROOT', Path.cwd() / 'CDRI_Odisha'))
    dirs = {k: root / v for k, v in {'raw': 'data/raw', 'external': 'data/external', 'interim': 'data/interim',
                                     'processed': 'data/processed', 'figures': 'outputs/figures',
                                     'tables': 'outputs/tables'}.items()}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    print('Project root:', root)
    return root, dirs


def union(gs):
    """Dissolve a GeoSeries/GeoDataFrame to one geometry (old and new geopandas)."""
    return gs.union_all() if hasattr(gs, 'union_all') else gs.unary_union


def need(path, notebook):
    """Stop with a clear message if an input from an earlier notebook is missing."""
    if not Path(path).exists():
        raise FileNotFoundError(f'{path} not found. Run Notebook {notebook} first.')
    return Path(path)


def gust_factor(exposure, avg_s):
    """WMO 3-s gust factor for each exposure class, given the averaging period (s) of the mean wind."""
    return np.array([GUST3[e][int(avg_s)] for e in np.asarray(exposure, dtype=object)], dtype=float)


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def _interp_valid(tn, t, y):
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() == 0:
        return np.full(len(tn), np.nan)
    if ok.sum() == 1:
        return np.full(len(tn), y[ok][0])
    return np.interp(tn, t[ok], y[ok])


def interp_track(t_h, lat, lon, vmax, rmw, pc=None, dt_h=1.0):
    """Interpolate a track to dt_h-hour steps and add the translation velocity (u east, v north; m/s)."""
    t_h = np.asarray(t_h, float)
    if len(t_h) < 2:
        return None
    tn = np.arange(t_h[0], t_h[-1] + 1e-9, dt_h)
    tr = {'t_h': tn, 'lat': _interp_valid(tn, t_h, lat), 'lon': _interp_valid(tn, t_h, lon),
          'vmax': _interp_valid(tn, t_h, vmax), 'rmw': _interp_valid(tn, t_h, rmw),
          'pc': _interp_valid(tn, t_h, pc) if pc is not None else np.full(len(tn), np.nan)}
    if len(tn) < 2:
        tr['u'] = tr['v'] = np.zeros(len(tn))
        return tr
    dy = np.gradient(tr['lat']) * 111.195
    dx = np.gradient(tr['lon']) * 111.195 * np.cos(np.radians(tr['lat']))
    dt = np.gradient(tn) * 3600.0
    tr['u'], tr['v'] = dx * 1000 / dt, dy * 1000 / dt
    return tr


def holland_max_wind(tr, plat, plon, p_env=1008.0, b_fixed=None, b_default=1.6, rho=1.15,
                     max_dist_km=600.0, chunk=4000):
    """Storm-maximum surface wind at points (same averaging basis as tr['vmax']).

    Holland (1980) radial profile scaled to the symmetric maximum wind, plus a translation asymmetry
    (Jelesnianski 1966: Vt·Rm·r/(Rm²+r²), half the translation speed at the RMW, strongest right of track).
    B comes from the central-pressure deficit (Holland 1980) when pressure is known, else b_default; clipped to 1.0–2.5.
    No land-roughness reduction is applied (as in STORM / Bloemendaal 2020); exposure is handled by the gust factor.
    """
    plat, plon = np.asarray(plat, float), np.asarray(plon, float)
    out = np.zeros(len(plat))
    lat, lon = tr['lat'][:, None], tr['lon'][:, None]
    u, v = tr['u'][:, None], tr['v'][:, None]
    vm, rm = tr['vmax'][:, None], np.maximum(tr['rmw'][:, None], 5.0)
    vt = np.hypot(u, v)
    vsym = np.maximum(vm - 0.5 * vt, 0.5 * vm)
    if b_fixed is not None:
        B = np.full_like(vm, float(b_fixed))
    else:
        dp = (p_env - tr['pc'][:, None]) * 100.0
        vg = vsym / SURF_TO_GRADIENT
        B = np.where(np.isfinite(dp) & (dp > 200.0), rho * np.e * vg ** 2 / np.maximum(dp, 200.0), b_default)
        B = np.clip(B, 1.0, 2.5)
    coslat = np.cos(np.radians(lat))
    for s in range(0, len(plat), chunk):
        pl, po = plat[None, s:s + chunk], plon[None, s:s + chunk]
        dy = (pl - lat) * 111.195
        dx = (po - lon) * 111.195 * coslat
        r = np.maximum(np.hypot(dx, dy), 0.5)
        x = (rm / r) ** B
        vs = vsym * np.sqrt(x * np.exp(1.0 - x))
        proj = (-dy * u + dx * v) / r                     # translation component along the cyclonic tangential direction
        w = vs + proj * (rm * r / (rm ** 2 + r ** 2))
        w = np.where(r <= max_dist_km, w, 0.0)
        out[s:s + chunk] = np.maximum(w.max(axis=0), 0.0)
    return out


def frag_p(im, theta, beta):
    """Lognormal fragility: P(damage | IM)."""
    im = np.maximum(np.asarray(im, float), 1e-3)
    return norm.cdf((np.log(im) - np.log(theta)) / beta)


def fit_fragility(im, y, theta0=None, beta0=0.2):
    """Maximum-likelihood lognormal fragility from binary outcomes (Baker 2015). Returns theta, beta, -loglik."""
    im, y = np.asarray(im, float), np.asarray(y, bool)
    ok = np.isfinite(im) & (im > 0)
    im, y = im[ok], y[ok]
    if theta0 is None:
        theta0 = np.median(im[y]) if y.any() else np.percentile(im, 99)
    li = np.log(im)

    def nll(p):
        pr = np.clip(norm.cdf((li - p[0]) / np.exp(p[1])), 1e-12, 1 - 1e-12)
        return -(y * np.log(pr) + (~y) * np.log(1 - pr)).sum()
    res = minimize(nll, [np.log(theta0), np.log(beta0)], method='Nelder-Mead',
                   options={'xatol': 1e-7, 'fatol': 1e-7, 'maxiter': 5000})
    return float(np.exp(res.x[0])), float(np.exp(res.x[1])), float(res.fun)


def bootstrap_fragility(im, y, n=200, seed=1):
    rng = np.random.default_rng(seed)
    im, y = np.asarray(im, float), np.asarray(y, bool)
    th, be = [], []
    for _ in range(n):
        i = rng.integers(0, len(im), len(im))
        if y[i].sum() < 3:
            continue
        t, b, _ = fit_fragility(im[i], y[i])
        th.append(t); be.append(b)
    return np.array(th), np.array(be)


# ---------- project data loaders ----------
def load_towers(dirs):
    """Towers from Notebook 00, with exposure from Notebook 01 if available (else a distance-to-coast rule)."""
    t = pd.read_parquet(need(dirs['processed'] / 'towers.parquet', '00'))
    ex = dirs['processed'] / 'towers_exposure.parquet'
    if ex.exists():
        e = pd.read_parquet(ex)
        t = t.merge(e[['tower_id', 'z0', 'terrain_cat', 'exposure', 'exposure_src']], on='tower_id', how='left')
    else:
        t['z0'], t['terrain_cat'], t['exposure'], t['exposure_src'] = np.nan, np.nan, None, None
    miss = t.exposure.isna()
    t.loc[miss, 'exposure'] = np.where(t.loc[miss, 'dist_coast_km'] <= 2.0, 'off-sea', 'in-land')
    t.loc[miss, 'exposure_src'] = 'rule: ≤2 km from coast = off-sea, else in-land'
    t['terrain_cat'] = t.terrain_cat.fillna(2).astype(int)
    return t


def load_network(dirs):
    """Lines and substations (with population served) from Notebooks 00/01."""
    import geopandas as gpd
    gpkg = need(dirs['processed'] / 'odisha_network.gpkg', '00')
    lines = gpd.read_file(gpkg, layer='lines')
    subs = gpd.read_file(gpkg, layer='substations')
    popf = dirs['processed'] / 'substations_population.csv'
    if popf.exists():
        p = pd.read_csv(popf)
        subs = subs.merge(p[['sub_id', 'pop_worldpop']], on='sub_id', how='left')
        subs['pop'] = subs.pop_worldpop.fillna(subs.pop_weight)
        subs['pop_src'] = np.where(subs.pop_worldpop.notna(), 'WorldPop 2020 service area', 'Raj et al. weight')
    else:
        subs['pop'], subs['pop_src'] = subs.pop_weight, 'Raj et al. weight'
    subs['pop'] = subs['pop'].fillna(0.0)
    return lines, subs


def load_fragility(dirs, which='any_damage'):
    f = dirs['processed'] / 'fragility_params.json'
    if f.exists():
        d = json.loads(f.read_text())
        if which in d:
            return d[which]['theta'], d[which]['beta'], d[which].get('source', 'Notebook 02')
    print('  ! fragility_params.json not found: using Raj et al. collapse fragility (292.5 km/h, β = 0.104)')
    return 292.5 / 3.6, 0.104, 'Raj et al. (collapse)'


# ---------- risk engine (Notebooks 05 and 06) ----------
def line_matrix(line_pos, n_lines):
    """Sparse tower × line indicator used to combine tower failure probabilities into line failure probabilities."""
    n = len(line_pos)
    return csr_matrix((np.ones(n), (np.arange(n), line_pos)), shape=(n, n_lines))


def storm_probabilities(H, cell, gfac, theta_t, beta, M, scale=1.0, batch=512):
    """For every storm: expected damaged towers, line failure probabilities and per-tower annual sums.

    H: storms × cells storm-max 10-min wind; cell: grid cell of each tower; gfac: wind→IM factor per tower;
    theta_t: fragility median per tower (array); M: tower × line indicator."""
    n_s, n_l = H.shape[0], M.shape[1]
    E = np.zeros(n_s); P_line = np.zeros((n_s, n_l), np.float32); tower_sum = np.zeros(len(cell))
    lt = np.log(theta_t)
    for s in range(0, n_s, batch):
        im = H[s:s + batch][:, cell].astype(np.float32) * (gfac * scale)[None, :]
        p = norm.cdf((np.log(np.maximum(im, 1e-3)) - lt[None, :]) / beta)
        E[s:s + batch] = p.sum(1)
        tower_sum += p.sum(0)
        logq = (M.T @ np.log1p(-np.minimum(p, 1 - 1e-9)).T).T
        P_line[s:s + batch] = 1 - np.exp(logq)
    return E, P_line, tower_sum


def outage_mc(P_line, frm, to, n_sub, sub_pop, n_mc=30, min_p=1e-3, seed=7):
    """Expected population disconnected per storm (Raj et al.: population outside the largest connected component)."""
    rng = np.random.default_rng(seed)

    def lost(alive):
        g = csr_matrix((np.ones(alive.sum()), (frm[alive], to[alive])), shape=(n_sub, n_sub))
        nc, lab = connected_components(g, directed=False)
        big = np.bincount(lab).argmax()
        return sub_pop[lab != big].sum()
    base = lost(np.ones(len(frm), bool))
    out = np.zeros(P_line.shape[0])
    for s in np.where(P_line.max(1) > min_p)[0]:
        fails = rng.random((n_mc, P_line.shape[1])) < P_line[s][None, :]
        out[s] = np.mean([max(lost(~f) - base, 0.0) for f in fails])
    return out


def annual_from_events(values, sim_year, n_years):
    return np.bincount(np.asarray(sim_year, int), weights=np.asarray(values, float), minlength=n_years)


def rp_table(annual, rps=(10, 25, 50, 100, 250, 500)):
    return {T: float(np.quantile(annual, 1 - 1 / T)) for T in rps if T <= len(annual)}


def exceedance_curve(annual):
    a = np.sort(np.asarray(annual))[::-1]
    rank = np.arange(1, len(a) + 1)
    return pd.DataFrame({'annual_value': a, 'return_period_yr': len(a) / rank})
def build_risk_context(dirs, storm_10min_corr=1.0, kv_scope=None):
    """Load everything the risk engine needs: STORM hazard (03), towers with exposure (00/01), lines and population."""
    z = np.load(need(dirs['interim'] / 'hazard_storm_cells.npz', '03'), allow_pickle=False)
    th = pd.read_parquet(need(dirs['processed'] / 'towers_hazard.parquet', '03'))[['tower_id', 'cell']]
    towers = load_towers(dirs).merge(th, on='tower_id', how='inner')
    if kv_scope:
        towers = towers[towers.voltage_kv.isin(kv_scope)]
    lines, subs = load_network(dirs)
    lines = lines.sort_values('line_id').reset_index(drop=True)
    towers['line_pos'] = towers.line_id.map({lid: i for i, lid in enumerate(lines.line_id)})
    towers = towers.dropna(subset=['line_pos']).reset_index(drop=True)
    towers['line_pos'] = towers.line_pos.astype(int)
    n_sub = int(max(subs.sub_id.max(), lines.from_sub.max(), lines.to_sub.max()))
    sub_pop = np.zeros(n_sub)
    sub_pop[subs.sub_id.values.astype(int) - 1] = subs['pop'].values
    return dict(H=z['H'], years=z['sim_year'].astype(int), ids=z['storm_id'], n_years=int(z['n_years']),
                towers=towers, lines=lines, subs=subs, cell=towers.cell.values.astype(int),
                gfac=gust_factor(towers.exposure, 600) * storm_10min_corr, M=line_matrix(towers.line_pos.values, len(lines)),
                frm=lines.from_sub.values.astype(int) - 1, to=lines.to_sub.values.astype(int) - 1, n_sub=n_sub, sub_pop=sub_pop)


def risk_run(ctx, theta_t, beta, scale=1.0, n_mc=30, cost_cr=0.73, outage_days=7.0, kwh_pp_day=3.0, voll_inr=100.0, seed=7):
    """One full risk run over all STORM storms. Costs in INR crore. Returns events, annual totals and EAD."""
    theta_t = np.broadcast_to(np.asarray(theta_t, float), (len(ctx['cell']),))
    E, P_line, tsum = storm_probabilities(ctx['H'], ctx['cell'], ctx['gfac'], theta_t, beta, ctx['M'], scale)
    pop = outage_mc(P_line, ctx['frm'], ctx['to'], ctx['n_sub'], ctx['sub_pop'], n_mc=n_mc, seed=seed)
    ev = pd.DataFrame({'storm_id': ctx['ids'], 'sim_year': ctx['years'], 'exp_damaged_towers': E, 'exp_pop_disconnected': pop})
    ev['tower_cost_cr'] = E * cost_cr
    ev['ens_mwh'] = pop * outage_days * kwh_pp_day / 1000.0
    ev['outage_cost_cr'] = ev.ens_mwh * 1000.0 * voll_inr / 1e7
    ev['total_cost_cr'] = ev.tower_cost_cr + ev.outage_cost_cr
    cols = ['exp_damaged_towers', 'exp_pop_disconnected', 'tower_cost_cr', 'outage_cost_cr', 'total_cost_cr']
    annual = pd.DataFrame({c: annual_from_events(ev[c], ev.sim_year, ctx['n_years']) for c in cols})
    return dict(events=ev, annual=annual, P_line=P_line, tower_rate=tsum / ctx['n_years'],
                ead={c: float(annual[c].mean()) for c in cols})
# ==== end of shared functions ====
