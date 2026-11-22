# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause

from pathlib import Path

project = "potlibc"
author = "Rohit Goswami"
copyright = "2026, Rohit Goswami"
version = "0.2.2"
release = version

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx_sitemap",
]

breathe_projects = {
    "potlibc": str(Path(__file__).parent / "_doxygen" / "xml"),
}
breathe_default_project = "potlibc"
breathe_default_members = ("members", "undoc-members")
breathe_domain_by_extension = {"h": "c"}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "rgpot": ("https://rgpot.rgoswami.me", None),
}

html_theme = "shibuya"
html_title = "potlibc"
html_baseurl = "https://omnipotentrpc.github.io/potlibc/"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "github_url": "https://github.com/OmniPotentRPC/potlibc",
    "accent_color": "teal",
    "dark_code": True,
    "globaltoc_expand_depth": 1,
    "nav_links": [
        {
            "title": "Ecosystem",
            "children": [
                {
                    "title": "rgpot",
                    "url": "https://rgpot.rgoswami.me",
                    "summary": "Potential-energy library and RPC server",
                },
                {
                    "title": "potlibc",
                    "url": "https://omnipotentrpc.github.io/potlibc/",
                    "summary": "C ABI for the PotLib 2023 surface catalog",
                },
            ],
        }
    ],
}

html_context = {
    "source_type": "github",
    "source_user": "OmniPotentRPC",
    "source_repo": "potlibc",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}
html_sidebars = {
    "**": [
        "sidebars/localtoc.html",
        "sidebars/repo-stats.html",
        "sidebars/edit-this-page.html",
    ]
}
sitemap_url_scheme = "{link}"

nitpicky = True
nitpick_ignore = [
    ("c:identifier", "size_t"),
    ("c:identifier", "uint32_t"),
    ("c:identifier", "int64_t"),
]
exclude_patterns = ["_doxygen"]

# These publisher domains return 403 to scripted HTTP clients (confirmed
# with both a default and a browser User-Agent), not just to Sphinx, so
# linkcheck cannot verify them. The DOIs that redirect here are instead
# checked against the Handle System API by tools/check_dois.py.
linkcheck_ignore = [
    r"^https://pubs\.aip\.org/",
    r"^https://aip\.scitation\.org/",
    r"^https://iopscience\.iop\.org/",
    # These links point at vendored files under upstream/potlib2023/ in this
    # repository. tests/test_catalog.py checks the catalog against those files.
    r"^https://github\.com/OmniPotentRPC/potlibc/blob/main/upstream/potlib2023/",
]
linkcheck_timeout = 20
linkcheck_retries = 2
