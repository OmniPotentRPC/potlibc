# SPDX-FileCopyrightText: 2026 Rohit Goswami
# SPDX-License-Identifier: BSD-3-Clause

from pathlib import Path

project = "potlibc"
author = "Rohit Goswami"
copyright = "2026, Rohit Goswami"
version = "0.1.0"
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
]
exclude_patterns = ["_doxygen"]
