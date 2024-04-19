from setuptools import setup, find_packages

with open("magpylib_force/__init__.py") as handle:
    for line in handle:
        if "__version__" in line:
            version = line.split(" = ")[-1].replace('"', "")
            break

with open("./README.md") as handle:
    readme_text = handle.read()

_short_description = (
    "An extension to the Magpylib library, providing force computation"
)

with open("./requirements.txt") as handle:
    requirements = [lr.strip() for lr in handle.read().splitlines() if lr.strip()]

setup(
    name="magpylib-force",
    version=version,
    description=_short_description,
    long_description=readme_text,
    long_description_content_type="text/markdown",
    author="Alexandre Boisselet",
    author_email="magpylib@gmail.com",
    url=("https://github.com/" "magpylib/magpylib-material-response"),
    license="MIT",
    packages=find_packages(),
    # include anything specified in Manifest.in
    install_requires=requirements,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords="",
)
