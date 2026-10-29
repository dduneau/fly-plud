# Fly PLUD

**Version 1.0.** A web app for the pathogen load upon death (PLUD) of infected flies, as described in the protocol (section "Analysis Bacterial Load Upon Death"). It screens each group for flies that did not die of the infection (Rosner test), compares the groups with a linear mixed model of log10 PLUD, and gives the plot, a methods paragraph and the R code.

The analysis is R itself (packages `lme4` and `EnvStats`). The app is built with shinylive, so R runs in the browser: there is no server, and the data are never uploaded. The first visit downloads R and the packages, which takes a few seconds.

## The data

The **per-fly table saved by [CFU Drop Counter](https://dduneau.github.io/cfu-drop-counter/)** (Load step): one row per fly, with its log10 load, its status, and the columns describing the fly (genotype, sex, box, day of injection, time of death …). The loads are estimated there, once, so that both apps use the same estimate.

Any table with one row per fly and a log10 (or log2) load works as well, separated by commas, semicolons, tabs or spaces, with decimal points or decimal commas (as saved by a French Excel). Column names are free: you tell the app which column plays which role. When a `Status` column is given, only flies with the status `Detected` are analysed; the others are counted and listed in a warning.

### Example data set

`example_PLUD_RAL821.csv` has the design the protocol recommends: three replicate experiments with 12 flies per group in each, giving 36 per group. The flies are DGRP line RAL821 infected with *P. rettgeri* at 35 days of age, from Duneau et al. (BMC Biology), sampled at random from the published data: 12 females and 12 males from each of three replicates. Females and males of a replicate shared one box, so the sexes are compared within boxes.

| Column | Content |
|---|---|
| `ID` | One identifier per fly, unique and with letters. |
| `Line`, `Sex`, `Age_at_infection` | What the fly is. |
| `Replicate` | The replicate experiment (day of injection). |
| `Box_ID` | The box, unique across the experiment and with letters (`Box_B1_35d`). Here one box per replicate, so the box also accounts for the replicate. With several boxes per replicate, give the replicate as the day of injection too. |
| `log10_CFU_per_fly`, `se_log`, `Status` | The load, its standard error and the status, as written by CFU Drop Counter (which also writes plates, drops and detection limits). |

## Using the app

The sidebar opens on the two steps every analysis needs; the other sections are optional and stay closed.

1. **Data.** Load the per-fly table, or try the example.
2. **Comparison.**
   - **Groups to compare** (x axis), for example `Sex` or `Genotype`.
   - **Which group is the least susceptible?** The app asks for it rather than choosing: the group with the lowest mortality in the survival assay. It is the reference, so that PLUD and mortality read in one direction (see the protocol). Differences are given relative to it.
   - **Order on the axis**: type the groups in the order you want, separated by commas. Groups left out of the list are added at the end. By default the order is alphabetical, with a group whose name starts with wild, WT, control or ctrl placed first.
   - **Box** (random effect). Tick **Box labels restart in each group** when box 1 of one line is not box 1 of another; the app ticks it when box numbers repeat across groups. **Day of injection** (optional) is a second random effect, for replicates holding several boxes.

**The screen** is the protocol's and is not adjustable: the Rosner test within each group (each combination of panel, x group and colour), at most k = 2 values for every group, a 5% threshold, low values only (a fly dying of another cause dies below its PLUD). High values that are flagged are shown, not removed. Results are given with and without the screen.

Optional sections:

- **Second factor and panels.** A **second factor, by colour**, shown side by side within each x group (for example age within each sex). Choose what to **compare** (the x groups within each colour, or the colours within each x group) and how to **fit**:
  - *one model for all groups*, for example `Sex + Sex:Age` when comparing ages within each sex: the groups share the residual variance and the box effects, and the likelihood ratio test is against `Sex` alone;
  - *a separate model for each group*, one test per sex.

  **Panels**: a third factor, for example the experiment. Each panel is screened and analysed separately.
- **Leave out groups.** Leave out whole groups, such as a genotype or a sex, without editing the file. Single flies cannot be removed by hand: only the screen removes flies, as described in the protocol.
- **Columns and labels.** The load column (log10 or log2; log2 is converted to log10), the **detection status** and the fly identifier, all found automatically in a table from CFU Drop Counter. CFU Drop Counter marks a fly `Detected` when it has a load; flies with another status (no colony, or every drop too crowded to count) are left out and counted in a message. Axis labels.

The app refuses settings that cannot be analysed: the box, the day or the fly identifier chosen as a group to compare, or a column with a single value. It warns when a group has fewer than 25 flies (the test is then less reliable; pool experimental replicates of the same group where possible), when high values are flagged, when each box holds a single group (the comparison then rests on the number of boxes), and when a variance is estimated at zero.

Females are drawn in orange and males in blue, as in the survival app.

## What it gives

- **Plots**: log10 PLUD per group (violins, every fly, the group mean, the number of flies, and the flies removed by the screen as crosses), and the difference of each group to the reference with its 95% confidence interval. PDF, PNG or SVG.
- **Results**, with and without the screen: the mean log10 PLUD of each group with its 95% CI; the difference to the reference with its profile-likelihood 95% CI and what it amounts to in bacteria per fly; the likelihood ratio test (maximum likelihood fits); the within-group SD, the SD of the random effects, and the difference this design detects with 80% power (2.8 × SD × √(2/n)).
- **Screen**: for each group, the number of flies and of values flagged and removed, and the list of flagged flies.
- **Flies**: the table of analysed flies with their group, box, log10 load, and whether the screen removed them.
- **Methods text**: a paragraph with the filters, the screen (with the number removed in each group), the model and the results.
- **R code**: a script that reproduces the screen and the model from the Flies table.

The screen finds values that do not fit a Gaussian distribution of log10 PLUD within a group. That they come from flies that died of something else is an inference, which the methods paragraph states. Do not apply it to loads measured in live flies (SPPL), where a low load means a fly controlling the infection.

## Repository layout and updating

```
app_plud/   the app: app.R and the example data set (edit here)
docs/       the published site, made by shinylive from app_plud/ (do not edit)
```

To update the app, edit `app_plud/app.R`, then in R, from the repository folder:

```r
shinylive::export("app_plud", "docs")
file.create("docs/.nojekyll")
httpuv::runStaticServer("docs")   # check it locally
```

The empty file `docs/.nojekyll` lets GitHub Pages serve all of shinylive's files. Then commit and push. GitHub Pages serves `docs/` from the `main` branch.

To run the app in R instead: `shiny::runApp("app_plud")`.

## Citation

If you use this app, please cite the protocol: [authors, title, journal, year, DOI].
