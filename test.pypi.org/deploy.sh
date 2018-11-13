# generate your package
python setup.py sdist
# upload your package to test pypi repository
twine upload --repository-url https://test.pypi.org/legacy/ dist/*

# intall your package form test.pypi
# pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple  hifise
