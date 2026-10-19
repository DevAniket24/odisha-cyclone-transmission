# Data sources

No raw data is stored in this repository. Notebook 00 downloads or creates everything in `MyDrive/CDRI_Odisha/data/`. All sources are free.

| Dataset | Used for | Source | Licence / terms |
|---|---|---|---|
| Odisha transmission network (41,814 towers, 227 substations, 327 lines) and Fani tower-damage records (128 towers, OPTCL via RTI) | Exposure, Fani case study, fragility | Raj, Kumar & Bhatia, [GitHub](https://github.com/surendervraj/cyclone-prepardness-for-regional-power-transmission-network) | No licence file; cite the authors. Tower-level damage data are not republished here. |
| STORM synthetic tropical cyclone tracks, present climate (10,000 years) | Hazard (Notebook 03) | Bloemendaal et al. 2020, [4TU.ResearchData](https://data.4tu.nl/articles/_/12706085) | CC BY 4.0 |
| IBTrACS v04r01, North Indian basin | Fani and validation cyclones; historical check | [NOAA NCEI](https://www.ncei.noaa.gov/products/international-best-track-archive) | Public domain |
| Natural Earth 10 m coastline | Distance to coast, coastal bands | [naturalearthdata.com](https://www.naturalearthdata.com) | Public domain |
| GADM 4.1 district boundaries (India level 2) | Districts | [gadm.org](https://gadm.org) | Free for non-commercial use |
| Meteostat hourly station data (NOAA ISD) | Fani wind check | [meteostat.net](https://meteostat.net) | CC BY-NC 4.0 |
| ESA WorldCover 10 m (Earth Engine) | Roughness / exposure class (Notebook 01) | ESA, via Google Earth Engine | CC BY 4.0 |
| WorldPop 100 m population 2020 (Earth Engine) | Population served by substations | WorldPop, via Google Earth Engine | CC BY 4.0 |
| Design codes: IS 875 (Part 3):1987, IS 802 (Part 1/Sec 1):1995 | Code comparison | [Public.Resource.Org](https://law.resource.org/pub/in/bis/S03/is.802.1.1.1995.pdf) | Public resource |
| IS 875 (Part 3):2015, IS 802 (Part 1/Sec 1):2015 | Code comparison (values to fill) | BIS, free download after registration | BIS terms |
| Cost and damage figures (Fani DLNA, CDRI, Raj et al., World Bank Lifelines, Entergy) | Cost-benefit inputs | See `REFERENCES.md` | As published |
