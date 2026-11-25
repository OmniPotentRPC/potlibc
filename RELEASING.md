<!-- SPDX-FileCopyrightText: 2026 Rohit Goswami -->
<!-- SPDX-License-Identifier: BSD-3-Clause -->

# Releasing potlibc

Cut a release from `main` after the changes since the latest `v*` tag are
in the tree. There is no fragment directory and no towncrier configuration.
Write the changelog from `git log <previous-tag>..HEAD` and the pull
requests merged in that range.

## Version strings

These six strings are the release version and must match:

- `meson.build`: the `version` keyword of `project()`. Meson writes that
  value into the installed library and `potlibc.pc`.
- `docs/source/conf.py`: `version`. `release` is assigned from it.
- `docs/Doxyfile`: `PROJECT_NUMBER`.
- `CITATION.cff`: `version`.
- `codemeta.json`: `version`.
- `.zenodo.json`: `version`.

Two dates follow the release date in the new changelog heading:
`date-released` in `CITATION.cff` and `dateModified` in `codemeta.json`.
`datePublished` in `codemeta.json` stays at the first release, 0.1.0.

Leave `meson_version` unchanged. It is the minimum Meson required to
configure the build, not the potlibc version. Leave historical changelog
headings unchanged.

## Steps

1. Add plain-prose entries under `## [Unreleased]` in `CHANGELOG.md`.
   Keep the existing [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
   sections. Describe the behavior a user can observe.

2. Replace the `## [Unreleased]` heading with `## [X.Y.Z] - YYYY-MM-DD`,
   using the date of the release commit. Do not leave an empty
   `## [Unreleased]` heading. The next change adds that heading back.

3. Set the six version strings above to `X.Y.Z` and the two release dates
   to `YYYY-MM-DD`. `CITATION.cff` must validate against the CFF 1.2.0
   schema, and both JSON files must parse:

   ```sh
   python -m json.tool codemeta.json > /dev/null
   python -m json.tool .zenodo.json > /dev/null
   ```

4. Search the work tree for the previous version. Remaining hits belong
   only in the historical changelog section or in prose that names that
   earlier release.

5. Run the `build` job from `.github/workflows/docs.yml`, except the
   artifact upload and GitHub Pages deploy steps. On a Debian or Ubuntu
   host the native packages are `doxygen`, `gfortran`, `liblapack-dev`,
   `meson`, and `ninja-build`, and Python is 3.12.

   ```sh
   python -m pip install -r docs/requirements.txt
   python docs/build.py
   python -m sphinx -b linkcheck -q docs/source docs/_build/linkcheck
   python tools/check_dois.py data/surfaces.json
   meson setup build
   meson compile -C build
   meson test -C build --print-errorlogs
   ```

6. Commit with the message `chore(release): X.Y.Z`. Author and committer
   are Rohit Goswami <rgoswami@ieee.org>. The message has no trailers.

7. After that commit is on `main`, create the annotated tag. Do not tag
   before the release commit is merged.

   ```sh
   git -c user.name='Rohit Goswami' -c user.email='rgoswami@ieee.org' \
     tag -a vX.Y.Z <release-commit> -m 'potlibc X.Y.Z'
   ```

   The annotated message is `potlibc` followed by the version, which is
   the form used for the previous tag.
