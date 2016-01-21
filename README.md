# hIPPYlib / web

This repo contains the [MkDocs](http://mkdocs.org) source files of the hIPPYlib website.

## Prerequisites

* Install [MkDocs](http://mkdocs.org)
    sudo pip install mkdocs

* Install [bootstrap](http://getbootstrap.com/) and [bootswatch](https://bootswatch.com/) themes
    sudo pip install mkdocs-bootstrap
    sudo pip install mkdocs-bootswatch


## To make changes to the website

* Clone this repository
    git clone https://github.com/hippylib/hippylib.github.io.git
* Edit or add some `.md` files (you may also need to update the `mkdocs.yml` config); 
* Preview locally
    mkdocs serve
* Publish on GitHub
    mkdocs gh-deploy
