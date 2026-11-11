# Authoring guide for OR Notes

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

## Figures

- Draw figures with matplotlib in `{python}` cells; code is hidden from readers
  (`echo: false` site-wide) but lives in the page.
- Import shared helpers and colors: `from _lib.orplot import new_axes, shade_polygon, convex_hull, mark_point, COLORS`.
  Colors: `p1` blue / `p2` orange for pieces or disjuncts; `weak` red for weaker
  relaxations; `hull` green for the convex hull or strongest relaxation.
- Every figure cell gets `#| label: fig-...` and a `#| fig-cap:` that says what to notice.
- Code runs with the project root as working directory (`execute-dir: project`).
- Python dependencies go in `requirements.txt`.

## Writing style

- Start with the plain-language idea and a tiny example, then the math, then why it works.
- Put proofs in collapsed callouts. Use `$...$` / `$$...$$` for math.
- Keep pages self-contained; link to related topics and paper pages.

## Build

- `quarto preview` for live preview, `quarto render` to build into `_site/`.
- Pushing to `main` publishes automatically via `.github/workflows/publish.yml`.
