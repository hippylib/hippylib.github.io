[![Build Status](https://travis-ci.org/hippylib/web.svg?branch=master)](https://travis-ci.org/hippylib/web)

# hIPPYlib / web

This repo contains the [MkDocs](http://mkdocs.org) source files of the hIPPYlib website.

## Edit the website directly from GitHub

The markdown sources for the hIPPYlib webpages are located in the `src` subfolder. Simply edit online the file you would like to update and commit. TravisCI will automatically build the html pages with `mkdocs` and push them to https://github.com/hippylib/hippylib.github.io.

> **NOTE**: This workflow is the easiest, however it will not allow you to preview changes in the website.

## Edit and preview the website locally

This will require installing some software on your computer, but it has the advantage that you'll be able to see the final results before publishing it on the web.

### Prerequisites

* Install [MkDocs](http://mkdocs.org)   
```
#!sh
    sudo pip install mkdocs==0.17.2
```

* Install [bootstrap](http://getbootstrap.com/) and [bootswatch](https://bootswatch.com/) themes
```
#!sh
    sudo pip install mkdocs-bootstrap==0.2.0
    sudo pip install mkdocs-bootswatch==0.5.0
```
    
* Install [MathJax](https://www.mathjax.org/) support using [python-markdown-math](https://github.com/mitya57/python-markdown-math)
```
#!sh
    sudo pip install python-markdown-math==0.6
```

Or, if you have conda installed, simply type

```conda-env create --file mkdocs-env.txt```



### To make changes to the website

* Clone this repository

```
#!sh
    git clone git@github.com:hippylib/web.git
```

* Edit the `.md` files located in the `src` folder. Note, if you add a new `.md` file you'll also need to update the `mkdocs.yml` config); 

* Preview locally
```
#!sh
    mkdocs serve
```

* Commit your changes and push them on GitHub
```
#!sh
    git add <files I want to commit>
    git commit -m "Explain changes"
    git push origin master:master
```
TravisCI integration will automatically call `mkdocs build` and push the freshly built *html pages* to https://github.com/hippylib/hippylib.github.io
