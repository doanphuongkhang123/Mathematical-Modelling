# Chemostat Beamer Presentation

This directory contains the unified English presentation for the chemostat
food-chain project. The main deck has 22 slides, followed by four backup
slides.

## Build

From the repository root:

```bash
make -C presentation
```

Output:

```text
presentation/chemostat_food_chain_presentation.pdf
```

Clean LaTeX build products with:

```bash
make -C presentation clean
```

## Source Layout

```text
main.tex                              document entry point
metadata.tex                          team and course metadata
preamble.tex                          Beamer theme and shared commands
sections/01_problem_and_approach.tex  motivation and assumptions
sections/02_model_and_theory.tex      model, equilibria, and stability
sections/03_numerical_method.tex      architecture, algorithms, and tests
sections/04_simulation_results.tex    P0-P3 and oscillatory experiments
sections/05_sensitivity_and_robustness.tex
sections/06_conclusions.tex
sections/appendix.tex
```

## Suggested Timing

| Content | Slides | Time |
|---|---:|---:|
| Problem and approach | 1-3 | 3.0 min |
| Model and theory | 4-9 | 7.0 min |
| Numerical method | 10-12 | 4.0 min |
| Simulation results | 13-17 | 6.5 min |
| Sensitivity and robustness | 18-20 | 4.0 min |
| Conclusions | 21 | 2.0 min |
| Q&A | 22 | 5.0 min |

The prepared presentation content is approximately 26.5 minutes before Q&A.
The appendix is for questions and is not part of the planned speaking time.

Replace the placeholders in `metadata.tex` with the five names and student
IDs before submission.
