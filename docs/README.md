# dcameaipe docs site

A static documentation/discussion site for the **dcameaipe** R package —
sidebar nav, search, numbered chapters, styled after the `interflex` manual.
No build step, no Jekyll — plain HTML/CSS/JS, meant to live inside the
[dcameaipe](https://github.com/SerranoSerrat/dcameaipe) repo itself.

## What's in here

```
index.html            Home page (D-CAME vs. AIPE diagram)
get-started.html      Ch. 1 — install & dcame_aipe() usage
theory.html           Ch. 2 — the estimands, and why linear models miss them
application.html      Ch. 3 — DGP2 examples + the robotization application
discussion.html       Ch. 4 — non-linear moderators, controls, discrete treatments
changelog.html        Release notes
references.html       Paper, package, and the surrounding literature
css/style.css         Shared styles
js/main.js            Sidebar search, active-link highlight, copy buttons
code/                 Standalone scripts linked from the chapters
img/                  Figures reproduced from the article
```

**Two files were renamed** from the previous version of this site. Delete the
old ones after copying:

```bash
git rm docs/applications.html docs/implications.html
```

`applications.html` → `application.html` and `implications.html` →
`discussion.html`.

## Assets still to add

**`img/`** — the chapters reference these filenames, taken from the article:

| File | Article figure |
|---|---|
| `img/figure1.png` | Fig. 1 — temperatures × SES, with exposure histograms |
| `img/figure2.png` | Fig. 2 — Monte Carlo of the linear interaction coefficient |
| `img/figure3.png` | Fig. 3 — estimated D-CAME (predicted values + CAMEs) |
| `img/figure4.png` | Fig. 4 — estimated AIPE (predicted values + IPEs + estimands) |
| `img/figure7.png` | Fig. 7 — the robotization application |

PNG at roughly 1600px wide renders well at the content column width. Until
they're added, the pages show broken-image placeholders with working captions.

**`changelog.html`** — the v0.1.0 release date is the only remaining `TODO`
in the site.

## Math rendering

`theory.html`, `application.html` and `discussion.html` load
[KaTeX](https://katex.org/) from a CDN (three tags in `<head>`, no build step,
works on GitHub Pages). `$…$` is inline and `$$…$$` is display. If you'd rather
not depend on a CDN, drop the three tags and the formulas degrade to plain
text — or self-host KaTeX under `css/`.

## Code

`code/kernel-conditional-skewness.R` is the simulation behind §4.2: it compares
the true CAME against the `interflex` kernel estimate under three treatment
distributions that share the same first two moments and the same correlation
with the moderator, differing only in conditional shape. Requires `interflex`
and `MASS`.

## 1. Add this to the existing `dcameaipe` repo, as `/docs`

```bash
# from inside your local clone of SerranoSerrat/dcameaipe
cp -r /path/to/this/folder docs

git add docs
git commit -m "Add documentation site"
git push
```

## 2. Point GitHub Pages at it

In the `dcameaipe` repo: **Settings → Pages → Source → Deploy from a
branch → `main` / `/docs` → Save**.

Live at:

```
https://serranoserrat.github.io/dcameaipe/
```

## 3. Link it from the README

Add near the top of `dcameaipe`'s main `README.md`:

```markdown
📖 **[Full documentation](https://serranoserrat.github.io/dcameaipe/)**
```

## Notes

- The sidebar is duplicated per page (no templating); active-link
  highlighting and search are computed client-side in `js/main.js`. If you
  add a page, copy the `<nav class="sidebar">…</nav>` block from an
  existing one and add a new `<li>`.
- The citation is now stated identically in `references.html` and should be
  reconciled with the "Citation" section of the package `README.md`, which
  previously used a different title.
- If you'd rather have the function-reference pages auto-generated and kept
  in sync with your roxygen comments in `man/`, look at
  [`pkgdown`](https://pkgdown.r-lib.org/) — it also builds into `/docs` and
  can coexist with hand-written "articles" for narrative content like
  Theory/Application/Discussion.
