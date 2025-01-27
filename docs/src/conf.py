import os
import sys
import toml
import shutil
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

sys.path.append(os.path.join(ROOT, "python"))
sys.path.append(os.path.join(ROOT, "docs", "extensions"))

# -- Project information -----------------------------------------------------

project = "Rascaline"
copyright = "2021, Rascaline developers"
author = "Rascaline developers"


def load_version_from_cargo_toml():
    with open(os.path.join(ROOT, "rascaline", "Cargo.toml")) as fd:
        data = toml.load(fd)
    return data["package"]["version"]


# The full version, including alpha/beta/rc tags
release = load_version_from_cargo_toml()


def build_cargo_docs():
    subprocess.run(["cargo", "doc", "--package", "rascaline", "--no-deps"])
    output_dir = os.path.join(ROOT, "docs", "build", "html", "reference", "rust")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    shutil.copytree(
        os.path.join(ROOT, "target", "doc"),
        output_dir,
    )


def extract_json_schema():
    subprocess.run(["cargo", "run", "--package", "rascaline-json-schema"])


def build_doxygen_docs():
    # we need to run a build to make sure the header is up to date
    subprocess.run(["cargo", "build", "--package", "rascaline-c-api"])
    subprocess.run(["doxygen", "Doxyfile"], cwd=os.path.join(ROOT, "docs"))


extract_json_schema()
build_cargo_docs()
build_doxygen_docs()


def setup(app):
    app.add_css_file(os.path.join(ROOT, "docs", "static", "rascaline.css"))


# -- General configuration ---------------------------------------------------

needs_sphinx = "3.5.0"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "breathe",
    "sphinx_tabs.tabs",
    "rascaline_json_schema",
    "html_hidden",
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["Thumbs.db", ".DS_Store"]


autoclass_content = "both"

breathe_projects = {
    "rascaline": os.path.join(ROOT, "docs", "build", "doxygen", "xml"),
}
breathe_default_project = "rascaline"
breathe_domain_by_extension = {
    "h": "c",
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]
