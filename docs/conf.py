import importlib.metadata

project = "protein-quest"
version = release = importlib.metadata.version(project)

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "autoapi.extension",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_immaterial",
    "sphinxcontrib.mermaid",
    "sphinx.ext.githubpages",
    "myst_nb",
    "cyclopts.sphinx_ext",
]

templates_path = ["_templates"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

language = "en"

html_title = "Protein Quest"
html_theme = "sphinx_immaterial"
html_static_path = ["_static"]
html_theme_options = {
    "icon": {
        "repo": "fontawesome/brands/github",
        "edit": "material/file-edit-outline",
    },
    "site_url": "https://www.bonvinlab.org/protein-quest/",
    "repo_url": "https://github.com/haddocking/protein-quest/",
    "repo_name": "protein-quest",
    "edit_uri": "blob/main/docs",
    "globaltoc_collapse": True,
    "features": [
        "navigation.expand",
        "navigation.tabs",
        # "navigation.tabs.sticky",
        # "toc.integrate",
        "navigation.sections",
        # "navigation.instant",
        # "header.autohide",
        "navigation.top",
        "navigation.footer",
        # "navigation.tracking",
        # "search.highlight",
        "search.share",
        "search.suggest",
        "toc.follow",
        "toc.sticky",
        "content.tabs.link",
        "content.code.copy",
        "content.action.edit",
        "content.action.view",
        "content.tooltips",
        "announce.dismiss",
    ],
    "palette": [
        {
            "media": "(prefers-color-scheme)",
            "toggle": {
                "icon": "material/brightness-auto",
                "name": "Switch to light mode",
            },
        },
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "teal",
            "accent": "teal",
            "toggle": {
                "icon": "material/lightbulb",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "teal",
            "accent": "teal",
            "toggle": {
                "icon": "material/lightbulb-outline",
                "name": "Switch to system preference",
            },
        },
    ],
    "version_dropdown": False,
    "toc_title_is_page_title": True,
    "social": [
        {
            "icon": "fontawesome/brands/github",
            "link": "https://github.com/haddocking/protein-quest",
            "name": "Source on github.com",
        },
        {
            "icon": "fontawesome/brands/python",
            "link": "https://pypi.org/project/protein-quest/",
        },
    ],
}

nb_execution_mode = "off"
myst_enable_extensions = ["colon_fence"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "aiohttp": ("https://docs.aiohttp.org/en/stable/", None),
    "dask": ("https://docs.dask.org/en/stable/", None),
    "distributed": ("https://distributed.dask.org/en/latest/", None),
    "yarl": ("https://yarl.aio-libs.org/en/latest/", None),
    "cyclopts": ("https://cyclopts.readthedocs.io/en/latest/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

always_document_param_types = True

autoapi_dirs = ["../src/protein_quest"]

# TODOS
# - mermaid render, does work in combined.py but not in README.md
# - solve warnings during build
# - links in docstrings
# - link from api doc to source code on GH
# - hide private functions/classes from docs
# - search working
# - cli toc should show sub-sub commands aka `filter combined`
