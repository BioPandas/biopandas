# How to Contribute ![](img/logos/1r8j_120.png)

I would be very happy about any kind of contributions that help to improve and extend the functionality of biopandas.

<br>
<br>

### Quick contributor checklist

- [ ]  Open a new "issue" on GitHub to discuss the new feature / bug fix  
- [ ]  Create and checkout a new topic branch   
- [ ]  Implement new feature or apply bugfix  
- [ ]  Add appropriate unit test functions  
- [ ]  Run `nosetests -sv` and make sure that all unit tests pass  
- [ ]  Add a note about the change to the `./docs/sources/CHANGELOG.md` file  
- [ ]  Modify documentation in `./docs/examples/` and `./docs/sources/`  
- [ ]  Push the topic branch to the server and create a pull request

<br>


## Getting Started

- If you don't have a [GitHub](https://github.com) account yet, please create one to contribute to this project.
- Please submit a ticket for your issue to discuss the fix or new feature before too much time and effort is spent for the implementation.

![](img/contributing/new_issue.png)

- Fork the `biopandas` repository from the GitHub web interface.

![](img/contributing/fork.png)

- Clone the `biopandas` repository to your local machine
	- `git clone https://github.com/<your_username>/biopandas.git`

<br>
<br>

## Making Changes

- Please avoid working directly on the master branch but create a new feature branch:
	- `git branch <new_feature>`
	- `git checkout <new_feature>`
- When you make changes, please provide meaningful `commit` messages:
	- `git add <edited_files>`
	- `git commit -m '<some note>'`
- Make an entry in the `biopandas/docs/sources/CHANGELOG.md` file.
- Add tests to the `biopandas/biopandas/tests` directory.
- Run all tests (e.g., via `nosetests`  or `pytest`).
- If it is a new feature, it would be nice (but not necessary) if you could update the documentation.
- Push your changes to a topic branch:
	- `git push -u origin <new_feature>`
- Submit a `pull request` from your forked repository via the GitHub web interface.

<br>

![](img/contributing/pull_request.png)

<br>
<br>

## Notes for Developers

### Building the documentation

Please note that documents containing code examples are generated from IPython Notebook files located in `mlxtend/docs/sources/ipynb` and converted to markdown via

```bash
~/github/mlxtend/docs/examples$ nbconvert --to markdown <file.ipynb>
```

The markdown file should be placed into the documentation directory at `biopandas/docs/sources` to build the documentation via  [mkdocs](http://www.mkdocs.org).
If you are adding a new document, please also include it in the pages section in the `biopandas/docs/mkdocs.yml` file.

To ensure that the documentation is rendered correctly, you can view the documentation locally by executing `mkdocs serve` from the `biopandas/docs` directory.

For example,

```bash
~/github/biopandas/docs$ mkdocs serve
```

The HTML files of the documentation are then build by executing

```bash
~/github/biopandas/docs$ mkdocs build
```

Deploying the documentation:

```bash
~/github/biopandas/docs$ mkdocs gh-deploy --clean
```

```bash
~/github/biopandas/docs$ mkdocs build --clean
~/github/biopandas/docs$ python setup.py upload_docs --upload-dir=site
```
