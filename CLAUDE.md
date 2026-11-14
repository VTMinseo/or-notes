# Authoring guide for OR Free Education

This repo is a Quarto website of Operations Research notes. Pages are written mostly by
Claude and reviewed by the owner. Follow these conventions so new pages slot in cleanly.

## Structure: category → subcategory → topic → paper

```
<category>/index.qmd                         # big category (order:, description:)
<category>/<subcategory>/index.qmd           # smaller category
<category>/<subcategory>/<topic>.qmd         # one topic per file
<category>/<subcategory>/papers/<author-year>.qmd   # optional paper notes
```

- Folder and file names: lowercase-kebab-case. Paper files: `firstauthor-year.qmd`.
- The sidebar is automatic (`contents: auto`). Every folder needs an `index.qmd`.
- Category/subcategory `index.qmd` pages use Quarto listings to list their children —
  copy an existing one (e.g. `integer-programming/index.qmd`).
- Order within a folder is set by the `order:` front-matter field.

## New pages

- Start from `_templates/topic.qmd` or `_templates/paper.qmd`.
- Required front matter: `kind`, `title`, `description`, `order` (topics), `categories`.
  Papers also need `author` and `year` (used in the papers table).
- `kind` is one of `category | subcategory | topic | paper`. Listings filter on it
  (Quarto globs are recursive), so a page without `kind` won't show up in any listing.
- Citations: add BibTeX entries to `references.bib`, cite with `[@key]`.
  Only add a DOI if it has been verified.

## Citation standard (write like an academic paper)

The site cites its sources the way a research paper would. These rules are not optional.

- **Credit the original paper.** Every named method, theorem, algorithm, cut family or
  model is attributed, at its first mention on each page, to the paper that introduced
  it — not only to a textbook or survey. Example: "Benders decomposition
  [@benders1962]", "lift-and-project cuts [@balas1993]", "Farkas' lemma [@farkas1902]".
  A survey may be cited *in addition*, never *instead*.
- **Name the authors in the text** when the paper is central to the page, using the
  in-text form `@key` ("@benders1962 showed …" renders as "Benders (1962) showed …").
- **Every reference must have a link.** Each `references.bib` entry needs a verified
  `doi` (preferred) or a stable `url` (arXiv abstract page, publisher page, report
  page). The reference list then shows authors, title, venue, year and a link.
- **Verify before adding.** Check authors, title, venue, volume, pages, year and DOI
  against the publisher, DOI registry, or the arXiv page. Never invent a DOI; if one
  can't be verified, use a verified URL instead.
- **Preprints** are cited as arXiv preprints with the arXiv ID, the version read, and
  the abstract-page URL (`https://arxiv.org/abs/<id>`), plus the arXiv DOI
  `10.48550/arXiv.<id>`.
- **Paper pages** open with the full citation box (authors, title, venue/arXiv ID,
  year, link) and are written in our own words: summarize, don't copy. Short quotes
  only, in quotation marks, with the source.
- **Figures adapted from a paper** say so in the caption ("in the spirit of Figure 3 of
  @fang2026"); figures are always redrawn from our own computations, never copied.

## Figures

- Draw figures with matplotlib in `{python}` cells; code is hidden from readers
  (`echo: false` site-wide) but lives in the page.
- Import shared helpers and colors: `from _lib.orplot import new_axes, shade_polygon, convex_hull, clip_halfplane, lattice_points, mark_point, COLORS`.
  `clip_halfplane` cuts a polygon along a split; `lattice_points` lists its integer points.
  Colors: `p1` blue / `p2` orange for pieces or disjuncts; `weak` red for weaker
  relaxations; `hull` green for the convex hull or strongest relaxation.
- Matrix-structure schematics: `from _lib.structure import draw_structure`, then
  `draw_structure(ax, "benders" | "dw", K=3)` (linking variables vs. linking constraints).
- Shared toy models, so related pages show consistent numbers: `_lib/toy_benders.py`
  (Benders loop, Routine 1 of Fang et al.) and `_lib/toy_dw.py` (column generation,
  LP / Dantzig–Wolfe / integer bounds).
- Every figure cell gets `#| label: fig-...` and a `#| fig-cap:` that says what to notice.
- Code runs with the project root as working directory (`execute-dir: project`).
- Python dependencies go in `requirements.txt`.

### Interactive figures

- Use Quarto's built-in Observable JS (`{ojs}` cells, with `d3` and `Inputs`); no
  server is needed. Compute the data in Python and pass it with
  `ojs_define(name=value)` (JSON-serializable lists/dicts only).
- The `ojs_define` cell must use `#| echo: false`, **not** `#| include: false`,
  otherwise the data never reaches the page.
- One `viewof` input per `{ojs}` cell. Use colors from `COLORS` (pass them through
  `ojs_define`) and `currentColor` for text, so dark mode works.
- Keep the static matplotlib version in a collapsed callout ("Static version of this
  figure (for printing)"). Example: the worked example in
  `integer-programming/disjunctive-programming/lift-and-project.qmd`.

## Writing style

- Start with the plain-language idea and a tiny example, then the math, then why it works.
- Put proofs in collapsed callouts. Use `$...$` / `$$...$$` for math.
- Keep pages self-contained; link to related topics and paper pages.

## Build

- `quarto preview` for live preview, `quarto render` to build into `_site/`.
- Pushing to `main` publishes automatically via `.github/workflows/publish.yml`.
