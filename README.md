# Fly PLUD

**Version 1.0.** A web app for the pathogen load upon death (PLUD) of infected flies, as described in the protocol (section "Analysis Bacterial Load Upon Death"). It screens each group for flies that did not die of the infection (Rosner test), compares the groups with a linear mixed model of log10 PLUD, and gives the plot, a methods paragraph and the R code.

The analysis is R itself (packages `lme4` and `EnvStats`). The app is built with shinylive, so R runs in the browser: there is no server, and the data are never uploaded. The first visit downloads R and the packages, which takes a few seconds.

## The data

The **per-fly table saved by [CFU Drop Counter](https://dduneau.github.io/cfu-drop-counter/)** (Load step): one row per fly, with its log10 load, its status, and the columns describing the fly (genotype, sex, box, day of injection, time of death …). The loads are estimated there, once, so that both apps use the same estimate.

Any table with one row per fly and a log10 (or log2) load works as well, separated by commas, semicolons, tabs or spaces. Column names are free: you tell the app which column plays which role. When a `Status` column is given, only flies with the status `Detected` are analysed; the others are counted and listed in a warning.

### Example data set

`example_PLUD_4lines_per_fly.csv`: 279 males of four lines (`Wild_type`, `Mutant_A`, `Mutant_B` and `Mutant_C`) collected at death after injection of *P. rettgeri*, two boxes per line, with the time of death. The table is in the format saved by CFU Drop Counter.

## Using the app

1. **Data.** Load the per-fly table or try the example. Choose the fly identifier, the load column (log10 or log2; log2 is converted to log10) and the status column.
2. **Flies to include.** **Keep only rows where** restricts the analysis to some values of a column (some genotypes, one sex). **A number lies between** restricts it to a range of a numeric column, for example flies that died at least 16 h after injection. Leave a bound empty for no limit.
3. **Design.**
   - **Groups on the x axis**, and their **order**: type the groups in the order you want, separated by commas (for example `Wild_type, Mutant_A, Mutant_B, Mutant_C`). Groups left out of the list are added at the end. By default the groups are in alphabetical order, with a group whose name starts with wild, WT, control or ctrl placed first.
   - **Second factor, by colour** (optional), shown side by side within each x group, for example a treatment. Then choose what to **compare**: the x groups within each colour, or the colours within each x group. And how to **fit**:
     - *one model for all groups*, for example `Line + Line:Treatment` when comparing treatments within lines: the groups share the residual variance and the box effects, and the likelihood ratio test is against `Line` alone, so it tests the treatment in all lines at once;
     - *a separate model for each group*, one test per line.
   - **Panels** (optional): a third factor, for example the experiment or the sex. Each panel is screened and analysed separately.
   - **Reference group**: differences are given relative to it. It is the first group of the order unless you choose another.
   - **Box** (random effect). Tick **Box labels restart in each group** when box 1 of one line is not box 1 of another; the app ticks it when box numbers repeat across groups. **Day of injection or experiment** (optional) is a second random effect.
4. **Screen.** The Rosner test is applied within each group (each combination of panel, x group and colour), with the same *k* for every group (2 by default) and a 5% threshold. By default only low values are removed: a fly dying of another cause dies below its PLUD. High values that are flagged are shown, not removed.

The app warns when a group has fewer than 25 flies (the test is then less reliable; pool experimental replicates of the same group where possible), when high values are flagged, when each box holds a single group (the comparison then rests on the number of boxes), and when a variance is estimated at zero.

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
