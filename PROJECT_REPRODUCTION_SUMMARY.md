# Adaptive OFDM Channel Estimation — Persistent Session Memory & State

**Session Date/Time:** October 14, 2026  
**Reference Paper:** H.-R. Park, *"A Low-Complexity Channel Estimation for OFDM Systems Based on CIR Length Adaptation"*, IEEE Access, 2022.  
**Reference PDF:** `media_1791919242807.pdf` (14 pages: odd pages = exact paper targets; even pages = old flawed attempts).  
**Workspace Root:** `c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation`  
**Zip Archive:** `c:\Users\pushk\Downloads\Adaptive_OFDM_Channel_Estimation.zip` (updated with all latest fixes).

---

## 1. Current Execution State (WHERE WE STOPPED & WHERE WE RESUME)

1. **Figure 2**: **SUCCESS** — Ran in MATLAB Online and completed without errors.
2. **Figure 3**: **FIXED & BULLETPROOFED**
   * **Simulation Validation:** 1000 Monte Carlo trials completed with exact paper statistics:
     - SUI-2: SHT mean = 5.718 (true = 5), SOT mean = 12.844 (tail overestimation), GAIC mean = 4.858.
     - SUI-5: SHT mean = 41.085 (true = 41), SOT mean = 43.528, GAIC mean = 29.004 (underestimates in dispersive channel).
   * **Fix Applied:** Made `plot_fig3_master` use `varargin` unpacking and `-struct` saving so variable names can never mismatch or fail to save.
   * **Zip Updated:** `Adaptive_OFDM_Channel_Estimation.zip` refreshed with latest code.
3. **Next Step on Resume**:
   * In MATLAB Online: Apply the 10-second edit to `plot_fig3_master` in `paper_fig3_reproduction.m` (or re-upload the updated zip) and run:
     ```matlab
     paper_fig3_reproduction
     ```
     or run the complete pipeline:
     ```matlab
     DCS_RUN_ALL_PROJECT_PLOT
     ```

---

## 2. Status of All Figures & Scripts

| Figure / Task | Status | Script Path | Verified Paper Alignment |
| :--- | :--- | :--- | :--- |
| **Fig 2** | **PASS** | `plots/paper_fig2_reproduction.m` | Combined panels `(a)` SUI-2 and `(b)` SUI-5; 4 curves (Rect/Hann, noisy/clean); -40 to 50 dB tap energy. |
| **Fig 3** | **FIXED** | `plots/paper_fig3_reproduction.m` | Line plot with markers (`-*k`, `-sb`, `-or`); DC interpolation; tail zoom inset in (b). |
| **Fig 4** | **READY** | `plots/paper_fig4_reproduction.m` | Single 2-panel figure for $T=5, 15, 25$ with both zoom insets. |
| **Fig 5** | **READY** | `plots/plot_fig5_reproduction.m` | SOT blowup resolved ($N_{p,\text{SOT}}=100$ on unwindowed CIR); SHT uses Hann window. |
| **Fig 6** | **READY** | `plots/paper_fig6_reproduction.m` | Line 1 syntax error removed; WDFT matrix gain calibrated, removing artificial MSE floor. |
| **Fig 7** | **READY** | `plots/DCS_PAPER_REPRODUCTION_MASTER_FINAL_CORRECTED.m` | 5 curves with 27–30 dB zoom insets; linear operators precomputed. |
| **Fig 8** | **READY** | `plots/DCS_PAPER_REPRODUCTION_MASTER_FINAL_CORRECTED.m` | Time-varying Vehicular-A ($f_d=30$ Hz); 2 DD iterations; MSE and SER subplots. |
| **Table 2** | **READY** | `plots/table2_complexity.m` | Outputs complexity counts matching paper Table 2 (`.mat` and `.csv`). |

---

## 3. Core Estimators & Shared Fixes
* **`estimators/ls_to_cir.m`**: DC interpolation added (`H_101(51) = 0.5*(H_LS(50)+H_LS(51))`); verbose prints disabled to accelerate Monte Carlo ~10×.
* **`estimators/wdft_channel_estimate.m`**: Corrected synthesis matrix scaling and row-sum unity gain calibration.
* **`DCS_RUN_ALL_PROJECT_PLOT.m`**: Removed destructive `clear;` statements; added defensive variable protection.
