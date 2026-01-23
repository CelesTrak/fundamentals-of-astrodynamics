# Instructions to generate the python documentation

To generate the documentation, install the necessary dependencies:

```
cd docs
pip install -r requirements.docs.txt
```

Then, you will need to register a Jupyter kernel so the notebooks can be included:

```
pip install ipykernel
python -m ipykernel install --user --name=vallado-env
```

Now you can generate the documentation using make

```
make html
```

You will find the HTML documentation under `_build/html/`. Open with `index.html` with your favourite browser!
