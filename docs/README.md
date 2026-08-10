# dcameaipe docs site

A static documentation/discussion site for the **dcameaipe** R package —
sidebar nav, search, numbered chapters, styled after the `interflex` manual.
No build step, no Jekyll — plain HTML/CSS/JS, meant to live inside the
[dcameaipe](https://github.com/SerranoSerrat/dcameaipe) repo itself.

## What's in here

```
index.html          Home page (dCAME vs. AIPE diagram)
get-started.html     Ch. 1 — install & real dcame_aipe() usage
applications.html    Ch. 2 — why the quantity of interest matters + compare="ALL"
implications.html    Ch. 3 — implications for research/review/policy
changelog.html        Release notes (skeleton — fill in as you tag releases)
references.html       Paper + package citation
css/style.css           Shared styles
js/main.js                Sidebar search, active-link highlight, copy buttons
```

Search for `TODO` across the files for what's still a placeholder
(mainly `changelog.html` and a couple of open-questions bullets in
`implications.html`).

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
- Consider reconciling the citation title used in this site's
  `references.html` / the package `README.md`'s "Citation" section — they
  currently differ slightly from the published article title.
- If you'd rather have the function-reference pages auto-generated and kept
  in sync with your roxygen comments in `man/`, look at
  [`pkgdown`](https://pkgdown.r-lib.org/) — it also builds into `/docs` and
  can coexist with hand-written "articles" for narrative content like
  Applications/Implications.
