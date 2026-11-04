# Pop Genome Mapper

Infer rough genomic positions of loci from population genotype data, using the linkage disequilibrium (LD) that neutral evolution leaves behind.

> **Status: early research prototype.** This README covers the *intended* tool. The code in the repo is a proof of concept for a simpler problem: mapping from a controlled testcross. See [Project status](#project-status) for what exists and what is planned. Sections marked *(planned)* describe interfaces that do not exist yet.

---

## Contents

- [The problem](#the-problem)
- [Why this is possible](#why-this-is-possible)
- [Method](#method)
- [Assumptions](#assumptions)
- [What the output means](#what-the-output-means)
- [Validation](#validation)
- [Installation](#installation)
- [Usage *(planned)*](#usage-planned)
- [Input format *(planned)*](#input-format-planned)
- [Repository layout](#repository-layout)
- [Project status](#project-status)
- [Roadmap](#roadmap)
- [Limitations](#limitations)
- [References](#references)
- [License](#license)

---

## The problem

**Input:** genotypes at $L$ biallelic loci for $n$ individuals sampled from one randomly mating population. Nothing about where the loci sit on the genome is known.

**Output:** for each locus, an estimate of:

1. which **linkage group** (chromosome) it belongs to, and
2. its rough **position** along that linkage group.

No pedigree, controlled cross, or reference genome is needed. The only signal is the statistical association between loci in the sampled population.

## Why this is possible

Recombination breaks up allele combinations every generation, and genetic drift in a finite population keeps creating new associations at random. Under neutrality the two forces reach a balance. Loci close together on a chromosome rarely recombine, so they stay correlated. Loci far apart, or on different chromosomes, become close to independent.

This correlation is **linkage disequilibrium**. For loci $A$ and $B$, with allele frequencies $p_A$, $p_B$ and haplotype frequency $p_{AB}$:

$$
D_{AB} = p_{AB} - p_A p_B, \qquad
r^2_{AB} = \frac{D_{AB}^2}{p_A(1-p_A)\,p_B(1-p_B)}
$$

At drift–recombination equilibrium, the expected $r^2$ falls steadily as the recombination fraction $c$ between the loci grows (Sved 1971):

$$
\mathbb{E}[r^2] \approx \frac{1}{1 + 4N_e c}
$$

$N_e$ is the effective population size. The Hill & Weir (1988) form also accounts for a finite sample of $n$ individuals:

$$
\mathbb{E}[r^2] = \frac{10 + C}{(2 + C)(11 + C)}\left[1 + \frac{(3 + C)(12 + 12C + C^2)}{n\,(2 + C)(11 + C)}\right], \qquad C = 4N_e c
$$

Because the curve only goes down, it can be inverted: observed LD gives an estimate of the scaled recombination rate $C$ for each pair of loci. Recombination rates behave like distances, and a set of pairwise distances can be turned into positions on a line.

> **Allele frequencies alone are not enough.** Single-locus allele frequencies carry no information about position. The tool needs **multi-locus genotypes per individual**, or at least two-locus genotype or haplotype frequencies, because the signal lives entirely in the *joint* distribution.

## Method

The pipeline has five stages. Each one is a separate, swappable module.

```
genotype matrix ──► 1. QC & filtering
                    2. pairwise LD (r²) ──► 3. linkage groups
                                           4. LD → distance
                                           5. 1-D embedding ──► positions per group
```

### 1. Quality control

- Drop monomorphic loci and loci below a minor-allele-frequency threshold. Rare alleles cap the attainable $r^2$ and make the estimates very noisy.
- Drop loci with too much missing data.
- Optionally test for Hardy–Weinberg equilibrium. Strong departures point to genotyping error or population structure.

### 2. Pairwise LD

Real data is usually **unphased** genotypes (0/1/2 copies of the alternate allele), not haplotypes. Haplotype $r^2$ is approximated by the squared Pearson correlation of genotype dosages (Rogers & Huff 2009). Phased haplotypes, when available, give $r^2$ directly.

The result is an $L \times L$ matrix of $\hat r^2$ values.

### 3. Linkage groups

For unlinked loci, $\mathbb{E}[r^2] \approx 1/n$, and $n\hat r^2$ is approximately $\chi^2_1$ distributed under the null of no association. The steps are:

1. Run the pairwise test and correct for multiple testing.
2. Build a graph whose edges join significantly linked pairs.
3. Split the graph into communities (e.g. Leiden or spectral clustering), not just connected components, so one spurious edge cannot merge two chromosomes.

Each community is treated as one linkage group.

### 4. LD to distance

Within each group, solve the Hill–Weir equation for $\hat C_{ij}$ for every pair. This gives a population-scaled recombination distance. Two practical issues come up:

- **Saturation:** once $\hat r^2$ falls near the $1/n$ floor, $\hat C$ becomes unbounded and carries no information. Those pairs are treated as "far" (censored) rather than as exact distances.
- **Additivity:** under a mapping function such as Haldane, $d = -\tfrac{1}{2}\ln(1 - 2c)$, the derived distances add up along the chromosome. That additivity is what the embedding step relies on.

### 5. One-dimensional embedding

Find positions $x_1, \dots, x_k$ for the $k$ loci in a group that best reproduce the estimated distances:

$$
\min_{x} \sum_{i<j} w_{ij}\,\bigl|\,|x_i - x_j| - \hat d_{ij}\,\bigr|
$$

The weights $w_{ij}$ reflect how informative each pair is. Close pairs get high weight and saturated pairs get low or zero weight. This is a *unidimensional scaling* problem, and the objective is **not convex**: each ordering of the loci forms its own basin. The intended approach is:

1. **Order:** find an initial ordering with a seriation method, such as spectral ordering on the LD matrix.
2. **Space:** with the order fixed, the absolute values open up into a linear program, solved exactly.
3. **Refine:** improve the order with local search, re-solving the spacing each time.

The solution is defined only up to **translation and reflection**. The tool fixes these explicitly: the first locus is placed at 0, and the orientation is chosen by convention. It does not pin positions with arbitrary bounds.

## Assumptions

The theory in [Why this is possible](#why-this-is-possible) holds under the following conditions. Violating them biases the results, usually toward too much LD between distant loci.

| Assumption | If violated |
|---|---|
| **Neutrality.** No selection at or near the loci | Selective sweeps inflate local LD, so loci look closer than they are |
| **Random mating** (panmixia) | Inbreeding or assortative mating inflates LD across the genome |
| **No population structure or recent admixture** | Creates LD between *unlinked* loci, which can merge linkage groups |
| **Drift–recombination equilibrium** | Bottlenecks or recent expansion shift the LD curve and bias $N_e$ |
| **Constant recombination rate** along the chromosome | Hotspots turn genetic distance into a warped version of physical distance |
| **Biallelic loci, accurate genotypes** | Genotyping error reduces $r^2$, so loci look farther apart |

## What the output means

- **Units.** LD identifies only $C = 4N_e c$, the population-scaled recombination rate. Without a separate estimate of $N_e$, positions are **relative**: order and spacing are meaningful, but absolute centimorgans are not. If $N_e$ is supplied, or two or more loci with known positions are given as anchors, positions can be converted to cM.
- **Precision.** LD estimates have high variance, especially at intermediate and long range. Expect the **order** of loci to be recovered well locally and the **spacing** only roughly. That is the "rough" in "rough genomic positions".
- **Orientation.** Which end of a linkage group is "0" is arbitrary.
- **Group labels.** Linkage-group numbering is arbitrary and does not match chromosome numbering.

## Validation

The method is only as credible as its recovery of a known truth. Every change is checked against **simulated populations where true positions are known**:

1. Simulate a neutral population with recombination using a coalescent simulator (e.g. [`msprime`](https://tskit.dev/msprime/)) or a forward Wright–Fisher simulator. Record the true map.
2. Sample $n$ individuals' genotypes.
3. Run the pipeline using the genotypes only.
4. Score the recovered map against the truth:
   - **Grouping:** adjusted Rand index between inferred and true linkage groups.
   - **Order:** Kendall's $\tau$ / Spearman's $\rho$ between inferred and true positions within each group, after fixing orientation.
   - **Spacing:** error in the inferred positions after the best translation, reflection and scale alignment.

The fit objective is **not** reported as an accuracy measure. A low objective value does not mean the map is correct.

Benchmarks will sweep $n$, $N_e$, number of loci, minor-allele-frequency cutoff, and departures from the assumptions above.

## Installation

The project uses [uv](https://docs.astral.sh/uv/) for reproducible environments.

```bash
git clone https://github.com/timlee-bioinf/Pop_Genome_Mapper.git
cd Pop_Genome_Mapper
uv sync          # (planned) creates .venv from pyproject.toml + uv.lock
```

Until `pyproject.toml` exists, the prototype notebooks need: `numpy`, `pandas`, `scipy`, `networkx`, `matplotlib`, `ipython`, `jupyter`.

## Usage *(planned)*

The commands below describe the target interface. They do not exist yet.

```bash
# Simulate a neutral population with a known map (for testing)
genome-mapper simulate --chromosomes 3 --loci 200 --ne 1000 --samples 500 \
    --seed 407 --out sim/

# Infer positions from genotypes
genome-mapper map sim/genotypes.csv --maf 0.05 --out results/

# Score an inferred map against the truth
genome-mapper evaluate results/positions.csv sim/truth.csv
```

```python
import genome_mapper as gm

sim = gm.simulate(chromosomes=3, loci=200, ne=1000, samples=500, seed=407)
result = gm.map_positions(sim.genotypes, maf=0.05)
report = gm.evaluate(result, sim.truth)
```

Every random step takes an explicit `seed`, so runs are reproducible from the command line alone.

## Input format *(planned)*

The input is a genotype matrix with **individuals as rows and loci as columns**, with values as alternate-allele dosage:

| value | meaning |
|---|---|
| `0` | homozygous reference |
| `1` | heterozygous |
| `2` | homozygous alternate |
| empty / `NA` | missing |

```csv
individual,locus_001,locus_002,locus_003
ind_1,0,1,2
ind_2,1,1,0
ind_3,2,NA,1
```

VCF input (`.vcf` / `.vcf.gz`) is planned so the tool can read standard variant-calling output directly.

Output is one row per locus: `locus, linkage_group, position, position_se`.

## Repository layout

**Current** (prototype):

| File | Purpose |
|---|---|
| `genome_generator.py` | Builds a random genome: chromosomes with genes at random bp positions |
| `meiosis_simu.py` | Simulates meiosis (≤1 crossover per chromosome) and fertilisation |
| `functions.py` | Testcross recombination fractions, linkage grouping, `basinhopping`-based position fitting |
| `solver.ipynb` | End-to-end demo: genome → 2,000 testcross offspring → inferred positions |
| `recombination_visualization.ipynb` | Plots chromosomes before and after crossover |
| `plots.ipynb`, `experiment.ipynb` | Exploration of the shape (non-convexity) of the fitting objective |
| `plots/` | Figures generated by the notebooks |

**Target:**

```
src/genome_mapper/
  simulate.py    # neutral population simulation (msprime wrapper + simple forward WF)
  io.py          # genotype matrix / VCF readers, result writers
  qc.py          # MAF, missingness, HWE filters
  ld.py          # r² from genotypes or haplotypes
  linkage.py     # significance testing + community detection → linkage groups
  distance.py    # Hill–Weir inversion, mapping functions, censoring
  embed.py       # seriation + LP spacing + local-search refinement
  evaluate.py    # comparison against known truth
  plot.py        # LD heatmaps, inferred-vs-true maps
  cli.py
tests/
notebooks/       # thin demos importing the package; outputs stripped
docs/
```

## Project status

| Component | Status |
|---|---|
| Genome simulation (fixed genes, bp positions) | ✅ prototype |
| Testcross offspring simulation | ✅ prototype |
| Recombination fraction from testcross | ✅ prototype |
| Linkage grouping (fixed threshold) | ✅ prototype, fragile |
| Position fitting (`basinhopping`, L1 loss) | ✅ prototype |
| Neutral *population* simulation | ⬜ planned |
| LD estimation from unphased genotypes | ⬜ planned |
| LD → distance (Hill–Weir inversion) | ⬜ planned |
| Statistical linkage grouping | ⬜ planned |
| Seriation + exact spacing solver | ⬜ planned |
| Evaluation against ground truth | ⬜ planned |
| Packaging, CLI, tests, CI | ⬜ planned |

The prototype answers a related but easier question: inferring positions from a **testcross**, where recombinants can be counted directly in one generation of offspring. The population-level problem swaps that direct count for an LD-based estimate. The embedding stage and the overall framing (pairwise distances → 1D positions) carry over unchanged.

## Roadmap

1. **Foundations.** Package layout, `pyproject.toml` + `uv.lock`, pytest, CI. Port the prototype into the package behind regression tests so nothing silently changes.
2. **Ground truth first.** Neutral population simulator with a known map, plus `evaluate.py`. Nothing is tuned before it can be scored.
3. **LD pipeline.** Genotype $r^2$ → Hill–Weir distances → embedding on a single linkage group.
4. **Linkage groups.** Statistical grouping across chromosomes.
5. **Better solver.** Seriation + LP spacing to replace `basinhopping`. Uncertainty via bootstrap over individuals.
6. **Robustness.** Benchmarks under structure, bottlenecks, hotspots and genotyping error.
7. **Real data.** VCF input, worked example on a public dataset.

## Limitations

- The LD signal decays fast. Beyond roughly $C \gg n$, distant pairs within a chromosome look unlinked. Large chromosomes may come out as several fragments unless enough loci sit densely along them.
- Absolute distances need $N_e$ or anchor loci (see [What the output means](#what-the-output-means)).
- Population structure is the most serious confounder. Run the tool on one homogeneous population, or correct for structure first.
- This is a research tool for rough positions. It does not replace a linkage map from a designed cross or a reference assembly.

## References

- Sved, J. A. (1971). Linkage disequilibrium and homozygosity of chromosome segments in finite populations. *Theoretical Population Biology*, 2(2), 125–141.
- Hill, W. G., & Weir, B. S. (1988). Variances and covariances of squared linkage disequilibria in finite populations. *Theoretical Population Biology*, 33(1), 54–78.
- Rogers, A. R., & Huff, C. (2009). Linkage disequilibrium between loci with unknown phase. *Genetics*, 182(3), 839–844.
- Haldane, J. B. S. (1919). The combination of linkage values and the calculation of distances between the loci of linked factors. *Journal of Genetics*, 8, 299–309.
- Kelleher, J., Etheridge, A. M., & McVean, G. (2016). Efficient coalescent simulation and genealogical analysis for large sample sizes. *PLoS Computational Biology*, 12(5), e1004842. (`msprime`)

## License

MIT. See [LICENSE](LICENSE).
