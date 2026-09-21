# Publication-table provenance

Each publication table has a repository generator and a documented chain back
to deposited inputs. Generated output directories are ignored by Git.

| Table | Generator | Immediate inputs | Upstream reproduction |
|---|---|---|---|
| Table 3 | `table_03/generate_table_03.m` | hydraulic statistics, case lookup, sweep limits, accepted-map summary, checkpoints | `code/data/videos/reproduce_hydraulic_statistics.m` and `reproduce_video_derived_csvs.m` |
| Table B1 | `table_b1/generate_table_b1.m` | synthetic autocorrelation summary | `code/data/syn/reproduce_synthetic_summary.m` |
| Table C1 | `table_c1/generate_table_c1.m` | accepted-map amplitude, Figure 8 velocity summary, Figure 9 depth summary | the real-video rebuild plus the Figure 8 and Figure 9 builders documented in those folders |

## Field-level origins

Table 3 hydraulic quantities are calculated from the deposited video
timestamps, stage/discharge series, long-term discharge values and surveyed
cross-section. Seeding density is calculated independently from the solver
checkpoint candidate counts, the accepted WSG-domain mask and sweep limits.

Table B1 prescribed values are parsed from each synthetic checkpoint filename.
Estimated wavelength and amplitude come from the selected WSG map. The velocity
difference is calculated using `U_s = sqrt(g*lambda/(2*pi))`.

Table C1 combines accepted autocorrelation wavelength and robust WSG amplitude
with the velocity and explicitly selected depth-root summaries. Its component
builders retain transect-level audit products so each case median can be traced
to the contributing values.
